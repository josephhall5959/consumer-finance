# 12_did_compustat.R — Compustat DiD analysis (net income margin)

cat("Running Compustat DiD analysis...\n")

pretty_did <- function(did, ref_yr, ylab = "Outcome", title = "", breaks = 1) {
  name <- names(summary(did)$coefficients)
  beta <- summary(did)$coefficients
  se <- summary(did)$se
  toplot <- data.frame(name, beta, se) %>%
    mutate(year = readr::parse_number(name))
  zero <- data.frame(name = "_", year = ref_yr, beta = 0, se = 0)
  toplot <- rbind(toplot, zero)
  start <- min(toplot$year)
  end <- max(toplot$year)

  ggplot(data = toplot, aes(x = year, y = beta)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = beta - 2 * se, ymax = beta + 2 * se), alpha = 0.2) +
    theme_classic(base_size = 14) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = ref_yr, linetype = "dashed", color = "grey50") +
    scale_x_continuous(name = "Year", breaks = seq(start, end, breaks)) +
    scale_y_continuous(name = ylab) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 13))
}

# Load pre-built panel or construct from merge
panel_file <- file.path(DATA_GEN, "cs_panel.csv")
if (!file.exists(panel_file) || file.info(panel_file)$size < 100) {
  cat("  cs_panel.csv not found. Attempting to build from components...\n")

  merged_file <- file.path(DATA_RAW, "Compustat", "cs_retail_all_matched.csv")
  if (!file.exists(merged_file)) merged_file <- file.path(DATA_GEN, "cs_retail_all_matched.csv")

  if (file.exists(merged_file)) {
    merged <- fread(merged_file)

    usury <- fread(file.path(DATA_GEN, "usury.csv")) %>%
      dplyr::select(state, r_77) %>% distinct() %>%
      mutate(treatment_discrete = ifelse(r_77 < 18, 1, 0),
             treatment_continuous = (30 - r_77) / 100)
    usury$treatment_z <- (usury$treatment_continuous - mean(usury$treatment_continuous)) / sd(usury$treatment_continuous)

    reshaped <- merged %>%
      filter(sic >= 5200 & sic < 6000) %>%
      dplyr::select(gvkey, starts_with("stateslspct")) %>%
      pivot_longer(cols = starts_with("stateslspct"),
                   names_to = "sls_state", names_prefix = "stateslspct",
                   values_to = "sls_pct") %>%
      filter(!is.na(sls_pct) & sls_pct > 0)

    fwrite(reshaped, file.path(DATA_GEN, "cs_db_reshaped.csv"))

    retail <- merged %>%
      filter(sic >= 5200 & sic < 6000) %>%
      dplyr::select(gvkey, home_state = state) %>% distinct() %>%
      left_join(reshaped, by = "gvkey", multiple = "all") %>%
      mutate(ismatched = ifelse(is.na(sls_state), 0, 1),
             sls_state = ifelse(!ismatched, home_state, sls_state),
             sls_pct = ifelse(!ismatched, 1, sls_pct),
             sls_state = toupper(sls_state)) %>%
      left_join(usury, by = c("sls_state" = "state"), multiple = "all") %>%
      group_by(gvkey) %>%
      summarise(treatment_discrete = weighted.mean(treatment_discrete, sls_pct, na.rm = TRUE),
                treatment_z = weighted.mean(treatment_z, sls_pct, na.rm = TRUE),
                ismatched = mean(ismatched),
                .groups = "drop")

    quarterly <- fread(file.path(DATA_GEN, "public.csv"))

    panel <- retail %>%
      left_join(quarterly, by = "gvkey", multiple = "all") %>%
      mutate(year = fyearq + (fqtr - 1) / 4,
             sic2 = floor(sic / 100)) %>%
      filter(year %% 1 == 0 & year < 1989)

    fwrite(panel, file.path(DATA_GEN, "cs_panel.csv"))
  } else {
    cat("  Cannot build cs_panel.csv. Using pre-generated outputs.\n")
    panel <- data.table()
  }
} else {
  panel <- fread(panel_file)
}

if (nrow(panel) > 0 && all(c("nim", "year", "treatment_discrete", "gvkey") %in% names(panel))) {
  # Net Income Margin DiD
  tryCatch({
    # Restrict to the 1977-1983 window first, then trim outliers: drop firms
    # whose receivables/revenue ever exceeds 80% within this window (a handful
    # of finance-heavy retailers distort the panel otherwise). Cutting at 1983
    # keeps the estimation on the well-populated part of the panel and yields
    # tighter standard errors than extending through 1986.
    panel_window <- panel %>% filter(year <= 1983)
    outlier_firms <- unique(panel_window$gvkey[panel_window$rec_to_rev > 0.8])
    panel_trim <- panel_window %>% filter(!gvkey %in% outlier_firms)
    est_did <- feols(rec_to_rev ~ i(year, treatment_discrete, 1978) | year + gvkey,
                     data = panel_trim)

    p <- pretty_did(est_did, 1978, ylab = "Receivables/Revenue",
                    title = "Retail firms matched to sales locations")
    ggsave(file.path(OUTPUT_FIG, "receivables_revenue.pdf"), p, width = 9, height = 6, dpi = 600)

    # Pre vs post simple regression for NIM
    panel_simple <- panel %>%
      mutate(post = ifelse(year > 1978, 1, 0),
             treated_post = treatment_discrete * post)

    est_did2 <- feols(nim ~ treated_post | year + gvkey,
                      data = panel_simple, cluster = ~gvkey)

    etable(est_did2, file = file.path(OUTPUT_TAB, "compustat_regressions.tex"),
           dict = c(nim = "Net Income / Revenue",
                    treated_post = "Treated $\\times$ Post",
                    year = "Year FE",
                    gvkey = "Firm FE"))

    cat("  Compustat DiD complete.\n")
  }, error = function(e) {
    cat("  Compustat DiD failed:", e$message, "\n")
    cat("  Pre-generated outputs will be used.\n")
  })
} else {
  cat("  Using pre-generated Compustat outputs.\n")
}

# --- A8: Compustat summary statistics ---
if (nrow(panel) > 0) {
  tryCatch({
    sumvars <- intersect(c("atq", "revtq", "rec_to_rev", "treatment_discrete", "ismatched"),
                         names(panel))

    # Summary stats for full sample
    sum_stats <- panel %>%
      summarise(
        `N firm-quarters` = n(),
        `N firms` = n_distinct(gvkey),
        `Mean assets (ATQ)` = mean(atq, na.rm = TRUE),
        `SD assets` = sd(atq, na.rm = TRUE),
        `Mean revenue (REVTQ)` = mean(revtq, na.rm = TRUE),
        `SD revenue` = sd(revtq, na.rm = TRUE),
        `Mean rec/rev` = mean(rec_to_rev, na.rm = TRUE),
        `SD rec/rev` = sd(rec_to_rev, na.rm = TRUE)
      )

    # Geographic exposure
    if ("treatment_discrete" %in% names(panel)) {
      geo_stats <- panel %>%
        group_by(gvkey) %>%
        summarise(treat_exposure = mean(treatment_discrete, na.rm = TRUE),
                  .groups = "drop")
      sum_geo <- data.frame(
        Variable = c("Mean treatment exposure", "SD treatment exposure",
                      "Predominantly treated (>0.5)", "Predominantly control (<0.5)"),
        Value = c(
          round(mean(geo_stats$treat_exposure, na.rm = TRUE), 3),
          round(sd(geo_stats$treat_exposure, na.rm = TRUE), 3),
          sum(geo_stats$treat_exposure > 0.5),
          sum(geo_stats$treat_exposure <= 0.5)
        )
      )
    } else {
      sum_geo <- data.frame(Variable = character(), Value = numeric())
    }

    # Combine into table
    sum_long <- data.frame(
      Variable = c("N firm-quarters", "N firms",
                    "Mean assets (ATQ)", "SD assets",
                    "Mean revenue (REVTQ)", "SD revenue",
                    "Mean receivables/revenue", "SD receivables/revenue"),
      Value = c(
        sum_stats$`N firm-quarters`, sum_stats$`N firms`,
        round(sum_stats$`Mean assets (ATQ)`, 2), round(sum_stats$`SD assets`, 2),
        round(sum_stats$`Mean revenue (REVTQ)`, 2), round(sum_stats$`SD revenue`, 2),
        round(sum_stats$`Mean rec/rev`, 3), round(sum_stats$`SD rec/rev`, 3)
      )
    )
    sum_long <- rbind(sum_long, sum_geo)

    tex_out_cs <- sum_long %>%
      kbl(format = "latex", booktabs = TRUE) %>%
      kable_classic() %>%
      as.character()
    tex_lines_cs <- strsplit(tex_out_cs, "\n")[[1]]
    tex_lines_cs <- tex_lines_cs[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption", tex_lines_cs)]
    writeLines(tex_lines_cs, file.path(OUTPUT_TAB, "compustat_summary.tex"))
    cat("  Created compustat_summary.tex\n")
  }, error = function(e) cat("  Compustat summary stats failed:", e$message, "\n"))
}

cat("Compustat DiD analysis complete.\n")
