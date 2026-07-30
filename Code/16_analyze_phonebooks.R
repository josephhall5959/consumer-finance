# 16_analyze_phonebooks.R — Yellow Pages / D&B phonebook analysis
# Requires db_phonebooks_merged.csv (proprietary D&B data matched to phonebooks)

cat("Analyzing phonebook/D&B data...\n")

db_pb_file <- file.path(DATA_RAW, "DB", "db_phonebooks_merged.csv")
if (!file.exists(db_pb_file)) {
  stop("db_phonebooks_merged.csv not found at: ", db_pb_file)
}

db_phonebooks <- fread(db_pb_file)
cat(sprintf("  Loaded %d rows from db_phonebooks_merged.csv\n", nrow(db_phonebooks)))

format_pct <- function(x) paste0(round(100 * x, 0), "%")

# Load usury treatment
usury <- fread(file.path(DATA_GEN, "usury.csv"))
usury_treat <- usury %>%
  dplyr::select(state, r_77) %>% distinct() %>%
  mutate(treated = ifelse(r_77 < 18, 1, 0))

# ---- 1. Adoption analysis (in-phonebook firms) ----
in_books <- db_phonebooks %>%
  filter(covered_by_phonebooks != "0") %>%
  mutate(
    SIC1_2dig = substr(as.character(sic1), 1, 2),
    industry = case_when(
      SIC1_2dig == "52" ~ "Hardware",
      SIC1_2dig == "53" ~ "Department Stores",
      SIC1_2dig == "54" ~ "Food",
      SIC1_2dig == "55" ~ "Auto",
      SIC1_2dig == "56" ~ "Clothing",
      SIC1_2dig == "57" ~ "Furnishing",
      SIC1_2dig == "58" ~ "Restaurants",
      SIC1_2dig == "59" ~ "Misc."
    ),
    acceptance = ifelse(offers_credit != "0", 1, 0),
    single = ifelse(ultdunsno == "", 1, 0),
    log_sls = log(as.numeric(sls)),
    log_sls_imp = ifelse(is.na(log_sls), 0, log_sls)
  )

# Filter to competitive markets (at least one acceptor per sic1-city-year)
in_sub <- in_books %>%
  group_by(sic1, city, year) %>%
  summarise(acceptance = sum(acceptance), .groups = "drop") %>%
  filter(acceptance > 0) %>%
  dplyr::select(sic1, city, year) %>%
  left_join(in_books, by = c("sic1", "city", "year")) %>%
  mutate(fe = paste0(sic1, industry, city))

# Card acceptance by industry table (descriptive).
# Emit a bare tabular (no table/caption wrapper): main.tex inputs this inside
# its own table float and inside a \resizebox, so an inner table environment
# would trigger "Not in outer par mode".
# Rates are computed within COMPETITIVE markets (in_sub: sic1-city-year cells
# with >=1 acceptor) -- i.e. places where bank cards are actually present -- so
# the fractions reflect adoption where it is a live choice rather than being
# diluted by markets with no card presence. The catch-all "Misc." retail group
# is excluded, and rows are sorted by the number of accepting firms, which puts
# clothing (the dominant card-using sector) at the top.
ae_tex <- in_sub %>%
  filter(!is.na(industry) & industry != "Misc.") %>%
  group_by(industry) %>%
  summarise(`Fraction Accepting` = mean(acceptance),
            `Firms Accepting` = sum(acceptance), .groups = "drop") %>%
  arrange(desc(`Firms Accepting`)) %>%
  mutate(`Fraction Accepting` = paste0(format(round(100 * `Fraction Accepting`, 2), nsmall = 2), "%")) %>%
  rename(Industry = industry) %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic() %>%
  as.character()
ae_lines <- strsplit(ae_tex, "\n")[[1]]
ae_lines <- ae_lines[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption", ae_lines)]
writeLines(ae_lines, file.path(OUTPUT_TAB, "card_acceptance_by_industry.tex"))
cat("  Created card_acceptance_by_industry.tex\n")

# Adoption regression (logit with marginal effects)
lm1 <- tryCatch({
  fit <- feglm(acceptance ~ log_sls + single | fe,
               data = in_sub, cluster = ~industry,
               family = binomial(link = "logit"))
  ame <- tryCatch(
    marginaleffects::avg_slopes(fit, vcov = FALSE),
    error = function(e2) { cat("  AME note:", e2$message, "\n"); NULL }
  )
  if (!is.null(ame)) {
    cat("  Phonebook logit AME:\n")
    print(ame)
  }
  fit
}, error = function(e) {
  cat("  Logit failed, falling back to OLS:", e$message, "\n")
  feols(acceptance ~ log_sls + single | fe,
        data = in_sub, cluster = ~industry)
})

# Add age using YEARSTARTED from retail.csv
retail_db <- file.path(DATA_GEN, "retail.csv")
db_file <- file.path(DATA_GEN, "db_70-85.csv")

birth_years <- NULL
if (file.exists(retail_db) && file.info(retail_db)$size > 1000) {
  db <- fread(retail_db, select = c("DUNSNO", "YEARSTARTED"))
  birth_years <- db[, .(birth_year = min(as.numeric(YEARSTARTED), na.rm = TRUE)),
                     by = .(dunsno = as.numeric(DUNSNO))]
  birth_years <- birth_years[is.finite(birth_year)]
  rm(db); gc()
} else if (file.exists(db_file) && file.info(db_file)$size > 1000) {
  db <- fread(db_file, select = c("DUNSNO", "YEARSTARTED"))
  birth_years <- db[, .(birth_year = min(as.numeric(YEARSTARTED), na.rm = TRUE)),
                     by = .(dunsno = as.numeric(DUNSNO))]
  birth_years <- birth_years[is.finite(birth_year)]
  rm(db); gc()
}

if (!is.null(birth_years) && nrow(birth_years) > 0) {
  sub_ages <- in_sub %>%
    left_join(birth_years, by = "dunsno") %>%
    mutate(age = year - birth_year)
  lm2 <- tryCatch({
    feglm(acceptance ~ log_sls + single + age | fe,
          data = sub_ages, cluster = ~industry,
          family = binomial(link = "logit"))
  }, error = function(e) {
    cat("  Logit with age failed, falling back to OLS:", e$message, "\n")
    feols(acceptance ~ log_sls + single + age | fe,
          data = sub_ages, cluster = ~industry)
  })
  adoption_dict <- c(log_sls = "log(Sales)",
                      single = "Single-establishment",
                      age = "Firm age",
                      acceptance = "Card acceptance",
                      fe = "Industry $\\times$ City FE")
  tex_out <- capture.output(etable(lm1, lm2, tex = TRUE, dict = adoption_dict))
} else {
  adoption_dict <- c(log_sls = "log(Sales)",
                      single = "Single-establishment",
                      acceptance = "Card acceptance",
                      fe = "Industry $\\times$ City FE")
  tex_out <- capture.output(etable(lm1, tex = TRUE, dict = adoption_dict))
}
writeLines(tex_out, file.path(OUTPUT_TAB, "adoption_firms.tex"))
cat("  Created adoption_firms.tex\n")

# ---- 2. DiD: Sales growth by acceptance × treatment (ypdb_difdif) ----
# This analysis merges phonebook acceptance (1977) with D&B sales panel (1970-1985).
# The full D&B panel (retail.csv, 17M rows) must be merged with phonebook status,
# which was originally done on the Stanford Yen servers.
# We use the pre-generated output from that analysis.
ypdb_src <- file.path(DATA_RAW, "DB", "ypdb_difdif.tex")
if (file.exists(ypdb_src)) {
  # Copy the pre-generated regression table and relabel its column/row headers
  # in plain English. The log(sales+1) detail is moved to the table notes in
  # main.tex; order of substitutions matters (compound labels before substrings).
  ypdb <- readLines(ypdb_src, warn = FALSE)
  subs <- list(
    c("SLS\\_growth",                    "Sales growth"),
    c("& d &",                           "& Survival (years) &"),
    c("acceptance $\\times$ treated",    "Accepts card $\\times$ Treated"),
    c("acceptance",                      "Accepts bank card"),
    c("log(SLS0+1)",                     "Log sales, 1977"),
    c("treated",                         "Treated"),
    c("STATE",                           "State"),
    c("SIC1",                            "Industry (SIC)")
  )
  for (s in subs) ypdb <- gsub(s[[1]], s[[2]], ypdb, fixed = TRUE)
  writeLines(ypdb, file.path(OUTPUT_TAB, "ypdb_difdif.tex"))
  cat("  Created ypdb_difdif.tex (relabeled headers)\n")
} else {
  stop("ypdb_difdif.tex not found. This table requires merging phonebook acceptance ",
       "with the full D&B panel (retail.csv) — originally done on the Yen servers.")
}

# --- A11: Yellow Pages summary statistics ---
tryCatch({
  yp_summary_rows <- list()

  # Number of phonebooks by year and state
  if ("year" %in% names(db_phonebooks) && "state" %in% names(db_phonebooks)) {
    yp_by_year <- db_phonebooks %>%
      group_by(year) %>%
      summarise(n_books = n_distinct(city), n_firms = n(), .groups = "drop")
    yp_by_state <- db_phonebooks %>%
      group_by(state) %>%
      summarise(n_books = n_distinct(paste0(city, year)), n_firms = n(), .groups = "drop")
  } else if ("city" %in% names(db_phonebooks)) {
    yp_by_year <- db_phonebooks %>%
      mutate(year = as.integer(year)) %>%
      group_by(year) %>%
      summarise(n_books = n_distinct(city), n_firms = n(), .groups = "drop")
  }

  # Match rates
  total_firms <- nrow(db_phonebooks)
  in_books_n <- nrow(in_books)
  match_rate <- round(in_books_n / total_firms, 3)

  # Acceptance rate: computed within competitive markets (in_sub), i.e. places
  # where bank cards are actually present, rather than over all covered firm-
  # years (which dilutes the rate with markets that have no card presence).
  accept_rate <- mean(in_sub$acceptance, na.rm = TRUE)

  yp_stats <- data.frame(
    Variable = c("Total firm-year observations", "Firms in phonebook coverage",
                  "Match rate", "Bank card acceptance rate (competitive markets)",
                  "N unique cities", "N unique years"),
    Value = c(
      total_firms, in_books_n,
      round(match_rate, 3), round(accept_rate, 3),
      n_distinct(db_phonebooks$city),
      n_distinct(db_phonebooks$year)
    )
  )

  tex_out_yp <- yp_stats %>%
    kbl(format = "latex", booktabs = TRUE) %>%
    kable_classic() %>%
    as.character()
  tex_lines_yp <- strsplit(tex_out_yp, "\n")[[1]]
  tex_lines_yp <- tex_lines_yp[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption", tex_lines_yp)]
  writeLines(tex_lines_yp, file.path(OUTPUT_TAB, "yp_summary.tex"))
  cat("  Created yp_summary.tex\n")
}, error = function(e) cat("  YP summary stats failed:", e$message, "\n"))

# --- A11b: Representativeness — matched vs unmatched firms in covered cities ---
tryCatch({
  # Compare acceptors vs non-acceptors on observables (within competitive markets)
  repr <- in_books %>%
    filter(!is.na(log_sls) & is.finite(log_sls)) %>%
    group_by(acceptance) %>%
    summarise(
      `Mean log(Sales)` = round(mean(log_sls, na.rm = TRUE), 2),
      `Frac. single-estab.` = round(mean(single, na.rm = TRUE), 3),
      N = n(),
      .groups = "drop"
    ) %>%
    mutate(Group = ifelse(acceptance == 1, "Accepts cards", "Does not accept")) %>%
    dplyr::select(Group, `Mean log(Sales)`, `Frac. single-estab.`, N)

  # Compare phonebook-covered vs not-covered firms (full dataset)
  all_firms <- db_phonebooks %>%
    mutate(
      covered = ifelse(covered_by_phonebooks != "0", "In phonebook area", "Not in phonebook area"),
      log_sls = log(as.numeric(sls)),
      single = ifelse(ultdunsno == "", 1, 0)
    ) %>%
    filter(!is.na(log_sls) & is.finite(log_sls))
  repr_coverage <- all_firms %>%
    group_by(covered) %>%
    summarise(
      `Mean log(Sales)` = round(mean(log_sls, na.rm = TRUE), 2),
      `Frac. single-estab.` = round(mean(single, na.rm = TRUE), 3),
      N = n(),
      .groups = "drop"
    ) %>%
    rename(Group = covered)

  repr_all <- bind_rows(repr_coverage, repr)
  tex_repr <- repr_all %>%
    kbl(format = "latex", booktabs = TRUE) %>%
    kable_classic() %>%
    as.character()
  tex_repr_lines <- strsplit(tex_repr, "\n")[[1]]
  tex_repr_lines <- tex_repr_lines[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption", tex_repr_lines)]
  writeLines(tex_repr_lines, file.path(OUTPUT_TAB, "yp_representativeness.tex"))
  cat("  Created yp_representativeness.tex\n")
}, error = function(e) cat("  Representativeness table failed:", e$message, "\n"))

# --- A11d: County-level representativeness (Yellow Pages sample vs rest of US) ---
# Compares CBP county aggregates for counties covered by the USTD collection
# against all other US counties. Covered counties are identified via the
# YP->CBP location crosswalk (Raw/CBP/crosswalk.dta, keyed to the 1970 CBP
# vintage). We therefore use 1970 CBP, the year the crosswalk is built for.
tryCatch({
  xw_file  <- file.path(DATA_RAW, "CBP", "crosswalk.dta")
  cbp_file <- file.path(DATA_GEN, "cbp_sic2.csv")
  if (file.exists(xw_file) && file.exists(cbp_file)) {
    xw  <- haven::read_dta(xw_file)
    cov <- unique(data.table::as.data.table(xw)[, .(stcode, cocode)])
    covkey <- paste(cov$stcode, cov$cocode)

    cbp <- data.table::fread(cbp_file)
    cbp_yr <- cbp[year == 1970 & cocode > 0 & !is.na(sic2)]

    cty <- cbp_yr[, .(
      emp    = sum(emp, na.rm = TRUE),
      est    = sum(est, na.rm = TRUE),
      retail = sum(est[sic2 >= 52 & sic2 <= 59], na.rm = TRUE)
    ), by = .(stcode, cocode)]
    cty[, covered := paste(stcode, cocode) %in% covkey]

    agg <- cty[, .(
      Counties = .N,
      avg_emp  = mean(emp),
      avg_est  = mean(est),
      avg_ret  = mean(retail),
      tot_emp  = sum(emp),
      tot_est  = sum(est)
    ), by = covered]

    ins <- agg[covered == TRUE]
    out <- agg[covered == FALSE]
    fmt <- function(x) formatC(round(x), format = "d", big.mark = ",")

    cty_lines <- c(
      "\\begin{tabular}[t]{lrr}",
      "\\toprule",
      "Variable & In Sample & Not In Sample\\\\",
      "\\midrule",
      sprintf("Avg Employment ('000) & %s & %s \\\\",
              fmt(ins$avg_emp / 1e3), fmt(out$avg_emp / 1e3)),
      sprintf("Avg Establishments & %s & %s \\\\",
              fmt(ins$avg_est), fmt(out$avg_est)),
      sprintf("Avg Retail Establishments & %s & %s \\\\",
              fmt(ins$avg_ret), fmt(out$avg_ret)),
      sprintf("Total Employment (M) & %s & %s \\\\",
              formatC(ins$tot_emp / 1e6, format = "f", digits = 0),
              formatC(out$tot_emp / 1e6, format = "f", digits = 0)),
      sprintf("Total Establishments ('000) & %s & %s \\\\",
              fmt(ins$tot_est / 1e3), fmt(out$tot_est / 1e3)),
      sprintf("Counties & %s & %s \\\\", fmt(ins$Counties), fmt(out$Counties)),
      "\\bottomrule",
      "\\end{tabular}"
    )
    writeLines(cty_lines, file.path(OUTPUT_TAB, "yp_county_repr.tex"))
    cat("  Created yp_county_repr.tex\n")
  } else {
    cat("  Skipping yp_county_repr.tex (crosswalk.dta or cbp_sic2.csv missing)\n")
  }
}, error = function(e) cat("  County representativeness table failed:", e$message, "\n"))

# --- A11c: Acceptance rates over time (balanced panel of cities) ---
tryCatch({
  # Use balanced panel: cities present in all years to avoid composition effects
  city_year_counts <- in_books %>%
    group_by(city) %>%
    summarise(n_years = n_distinct(year), .groups = "drop")
  max_years <- max(city_year_counts$n_years)
  balanced_cities <- city_year_counts %>% filter(n_years == max_years) %>% pull(city)

  if (length(balanced_cities) >= 5) {
    accept_ts <- in_books %>%
      filter(city %in% balanced_cities) %>%
      group_by(year) %>%
      summarise(acceptance_rate = mean(acceptance, na.rm = TRUE),
                n_accepting = sum(acceptance, na.rm = TRUE),
                n_firms = n(), .groups = "drop") %>%
      mutate(year = as.integer(year))
  } else {
    # Fallback: all cities
    accept_ts <- in_books %>%
      group_by(year) %>%
      summarise(acceptance_rate = mean(acceptance, na.rm = TRUE),
                n_accepting = sum(acceptance, na.rm = TRUE),
                n_firms = n(), .groups = "drop") %>%
      mutate(year = as.integer(year))
  }

  p_ts <- ggplot(accept_ts, aes(x = year, y = acceptance_rate)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    theme_classic(base_size = 14) +
    scale_y_continuous(name = "Card acceptance rate", labels = scales::percent_format()) +
    scale_x_continuous(name = "Year") +
    theme(axis.title = element_text(size = 13))
  ggsave(file.path(OUTPUT_FIG, "acceptance_over_time.pdf"), p_ts, width = 7, height = 5)
  cat(sprintf("  Created acceptance_over_time.pdf (%d balanced cities)\n", length(balanced_cities)))
}, error = function(e) cat("  Acceptance time series failed:", e$message, "\n"))

rm(db_phonebooks, in_books, in_sub); gc()
cat("Phonebook analysis complete.\n")
