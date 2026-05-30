# 13_did_marquette.R — D&B firm-level DiD (Marquette decision)
# Uses db_state.csv (state-level aggregates) and retail.csv for recession analysis.

cat("Running Marquette DiD analysis...\n")

# Look for state-level aggregates first (memory-efficient)
state_file <- file.path(DATA_GEN, "db_state.csv")
retail_file <- file.path(DATA_GEN, "retail.csv")
if (!file.exists(retail_file)) retail_file <- file.path(DATA_RAW, "DB", "retail.csv")

# Load usury treatment assignment
usury <- fread(file.path(DATA_GEN, "usury.csv")) %>%
  dplyr::select(state = state, r_77) %>%
  distinct() %>%
  mutate(treated = ifelse(r_77 < 18, 1, 0))

# ---- 1. Event-study DiD using state-level aggregates ----
if (file.exists(state_file) && file.info(state_file)$size > 100) {
  cat("  Loading state-level aggregates...\n")
  state_agg <- fread(state_file)

  state_agg <- merge(state_agg, usury, by.x = "STATE", by.y = "state", all.x = TRUE)
  state_agg <- state_agg[!is.na(treated)]
  state_agg[, exit_rate := exit]

  cat("  Running exit rate event study (state-level)...\n")
  tryCatch({
    did_exit <- feols(exit_rate ~ i(YEAR, treated, 1977) | STATE + YEAR,
                      data = state_agg[YEAR >= 1970],
                      weights = ~n,
                      cluster = ~STATE)

    name <- names(coef(did_exit))
    beta <- coef(did_exit)
    se <- summary(did_exit)$se
    toplot <- data.frame(name, beta, se)
    toplot <- toplot[grep("treated", toplot$name), ] %>%
      mutate(year = readr::parse_number(name))
    zero <- data.frame(name = "_", year = 1977, beta = 0, se = 0)
    toplot <- rbind(toplot, zero)

    p <- ggplot(data = toplot, aes(x = year, y = beta)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      geom_ribbon(aes(ymin = beta - 2 * se, ymax = beta + 2 * se), alpha = 0.2) +
      theme_classic(base_size = 14) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 1977, linetype = "dashed", color = "grey50") +
      annotate("text", x = 1977.2, y = Inf, label = "Marquette", vjust = 2,
               color = "grey30", size = 3.5, hjust = 0) +
      scale_x_continuous(name = "Year", breaks = seq(1970, 1985, 1)) +
      scale_y_continuous(name = "Firm failure rate") +
      theme(axis.title = element_text(size = 13))
    ggsave(file.path(OUTPUT_FIG, "difdif_exit_all.pdf"), p, height = 6, width = 6, units = "in")
    cat("    Created difdif_exit_all.pdf\n")

    # A7: With state-level controls
    bea_ctrl_file <- file.path(DATA_GEN, "bea.csv")
    if (file.exists(bea_ctrl_file)) {
      tryCatch({
        bea_ctrl <- fread(bea_ctrl_file)
        bea_77_ctrl <- bea_ctrl %>%
          filter(year == 1977) %>%
          group_by(state) %>%
          summarise(baseline_wage = mean(wage_rate, na.rm = TRUE),
                    .groups = "drop")
        state_agg_ctrl <- merge(state_agg, bea_77_ctrl,
                                by.x = "STATE", by.y = "state", all.x = TRUE)

        did_exit_ctrl <- feols(exit_rate ~ i(YEAR, treated, 1977) |
                                 STATE + YEAR + YEAR[baseline_wage],
                               data = state_agg_ctrl[YEAR >= 1970],
                               weights = ~n, cluster = ~STATE)
        name_c <- names(coef(did_exit_ctrl))
        beta_c <- coef(did_exit_ctrl)
        se_c <- summary(did_exit_ctrl)$se
        toplot_c <- data.frame(name = name_c, beta = beta_c, se = se_c)
        toplot_c <- toplot_c[grep("treated", toplot_c$name), ] %>%
          mutate(year = readr::parse_number(name))
        zero_c <- data.frame(name = "_", year = 1977, beta = 0, se = 0)
        toplot_c <- rbind(toplot_c, zero_c)

        p_c <- ggplot(data = toplot_c, aes(x = year, y = beta)) +
          geom_line(linewidth = 0.8) + geom_point(size = 2) +
          geom_ribbon(aes(ymin = beta - 2 * se, ymax = beta + 2 * se), alpha = 0.2) +
          theme_classic(base_size = 14) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          geom_vline(xintercept = 1977, linetype = "dashed", color = "grey50") +
          scale_x_continuous(name = "Year", breaks = seq(1970, 1985, 1)) +
          scale_y_continuous(name = "Firm failure rate") +
          theme(axis.title = element_text(size = 13))
        ggsave(file.path(OUTPUT_FIG, "difdif_exit_controls.pdf"),
               p_c, height = 6, width = 6, units = "in")
        cat("    Created difdif_exit_controls.pdf\n")
        rm(bea_ctrl, bea_77_ctrl, state_agg_ctrl, did_exit_ctrl); gc()
      }, error = function(e) cat("    Controls regression failed:", e$message, "\n"))
    }

    rm(did_exit); gc()
  }, error = function(e) cat("    Exit DiD failed:", e$message, "\n"))

  rm(state_agg); gc()
} else {
  cat("  db_state.csv not available for event study.\n")
}

# ---- 1b. SMM target moments: sales DiD (feeds the structural model) ----
# Moment 1: average sales change from Marquette   (avg treated firm, post)
# Moment 2: accept-vs-not sales differential       (card interaction)
# We use the CONTINUOUS usury-tightness intensity gap = max(0, 18 - r_77) rather
# than a binary treated dummy: the gap exploits the full cross-state variation in
# how strongly Marquette relaxed each state's binding rate ceiling, which extracts
# more signal (lower moment SEs) than a 0/1 split.  Both moments come off a SINGLE
# interaction regression so we keep their full 2x2 joint covariance, then scale to
# the average binding gap among treated states so the moments retain the "effect
# for a typical treated state" interpretation that the structural model expects.
#   moment1 (avg effect)   = gbar * (b_gp + p_card * b_gp:card)
#   moment2 (differential) = gbar *  b_gp:card
# where gbar = mean(gap | treated).  Note the scaling cancels in the t-statistics,
# so the precision gain over the binary spec reflects genuine added signal.  The
# clustered estimates AND covariance are written to smm_moments.csv and read by
# 21_nested_logit.R, so structural-parameter uncertainty propagates from these
# reduced-form regressions through the SMM to SE(sigma)/SE(cost) via delta method.
if (file.exists(retail_file) && file.info(retail_file)$size > 1000) {
  cat("  Estimating SMM target moments from firm-level sales...\n")
  tryCatch({
    rs <- fread(retail_file,
                select = c("STATE", "YEAR", "SIC1", "DUNSNO",
                           "log_sls", "card", "treated", "r_77"))
    rs <- rs[!is.na(log_sls) & is.finite(log_sls) & !is.na(treated) & !is.na(r_77)]
    rs[, post := as.integer(YEAR >= 1978)]
    rs[, gap := pmax(0, 18 - r_77)]      # usury bindingness (0 for control states)
    rs[, gap_post := gap * post]

    m <- feols(log_sls ~ gap_post * card | DUNSNO + YEAR,
               data = rs, cluster = ~STATE)
    b <- coef(m); V <- vcov(m); nm <- names(b)
    i_g  <- which(nm == "gap_post")
    i_gc <- which(nm == "gap_post:card")

    # Average binding gap among treated states and card share among treated-post
    # firm-years -- the weights that turn per-unit-gap coefficients into the
    # average effect for a typical treated state.
    gbar   <- rs[treated == 1, mean(gap)]
    p_card <- rs[gap_post > 0, mean(card)]

    # Linear map L: moments = L %*% b  (scaled to gbar)
    L <- matrix(0, nrow = 2, ncol = length(b))
    L[1, i_g] <- gbar; L[1, i_gc] <- gbar * p_card
    L[2, i_gc] <- gbar
    moments <- as.numeric(L %*% b)
    Sigma_m <- L %*% V %*% t(L)          # full 2x2 covariance, off-diagonal kept

    smm <- data.frame(
      moment   = c("avg_sales_change", "accept_differential"),
      estimate = moments,
      se       = sqrt(diag(Sigma_m)),
      cov12    = Sigma_m[1, 2]           # covariance between the two moments
    )
    fwrite(smm, file.path(DATA_GEN, "smm_moments.csv"))
    cat(sprintf("    Avg treated usury gap (gbar):   %.2f pts\n", gbar))
    cat(sprintf("    Moment 1 (avg sales change):    %.4f (se %.4f)\n",
                moments[1], sqrt(Sigma_m[1, 1])))
    cat(sprintf("    Moment 2 (accept differential): %.4f (se %.4f)\n",
                moments[2], sqrt(Sigma_m[2, 2])))
    cat(sprintf("    Moment correlation:             %.3f\n",
                Sigma_m[1, 2] / sqrt(Sigma_m[1, 1] * Sigma_m[2, 2])))
    cat("    Created smm_moments.csv\n")
    rm(rs, m); gc()
  }, error = function(e) cat("    SMM moment estimation failed:", e$message, "\n"))
}

# ---- 2. Recession triple-difference ----
# Aggregate to state-SIC1-year-single cells to fit in memory, then run DiD.
if (file.exists(retail_file) && file.info(retail_file)$size > 1000) {
  cat("  Running recession triple-difference (cell-level aggregation)...\n")
  tryCatch({
    # Read only needed columns
    retail <- fread(retail_file,
                    select = c("STATE", "YEAR", "SIC1", "exit",
                               "log_sls", "single", "treated"))
    retail <- retail[!is.na(log_sls) & !is.na(exit)]
    retail[, SIC1 := as.integer(SIC1)]
    retail[, high_acceptance := fifelse(SIC1 >= 5500 & SIC1 < 5800, 1L, 0L, na = 0L)]

    # Aggregate to cells: state × SIC1 × year × single
    cells <- retail[, .(exit_rate = mean(exit, na.rm = TRUE),
                        mean_log_sls = mean(log_sls, na.rm = TRUE),
                        n = .N),
                    by = .(STATE, YEAR, SIC1, single, treated, high_acceptance)]
    rm(retail); gc()

    cat(sprintf("    Aggregated to %d cells from retail.csv\n", nrow(cells)))

    cells[, fe := paste0(STATE, SIC1, YEAR)]

    run_recession <- function(dt, rec_year) {
      dt[, recession := fifelse(YEAR == rec_year, 1L, 0L)]
      dt[, treated_recession := treated * recession]
      dt[, treated_single := treated * single]
      dt[, single_recession := single * recession]
      dt[, treated_recession_single := treated * recession * single]

      # Original composite FE
      lm_all <- feols(exit_rate ~ mean_log_sls + treated_recession + treated_single +
                        single_recession + treated_recession_single | fe,
                      data = dt, weights = ~n, cluster = ~STATE)
      lm_hi <- feols(exit_rate ~ mean_log_sls + treated_recession + treated_single +
                       single_recession + treated_recession_single | fe,
                     data = dt[high_acceptance == 1], weights = ~n, cluster = ~STATE)

      # A9: Proper three-way FE structure: STATE^YEAR + SIC1^YEAR + STATE^SIC1
      lm_all_3way <- tryCatch({
        feols(exit_rate ~ mean_log_sls + treated_recession + treated_single +
                single_recession + treated_recession_single |
                STATE^YEAR + SIC1^YEAR + STATE^SIC1,
              data = dt, weights = ~n, cluster = ~STATE)
      }, error = function(e) { cat("    3-way FE failed:", e$message, "\n"); NULL })

      lm_hi_3way <- tryCatch({
        feols(exit_rate ~ mean_log_sls + treated_recession + treated_single +
                single_recession + treated_recession_single |
                STATE^YEAR + SIC1^YEAR + STATE^SIC1,
              data = dt[high_acceptance == 1], weights = ~n, cluster = ~STATE)
      }, error = function(e) { cat("    3-way FE (hi) failed:", e$message, "\n"); NULL })

      list(lm_all, lm_hi, lm_all_3way, lm_hi_3way)
    }

    # 1983 recession
    res83 <- run_recession(cells, 1983)
    # 1975 recession placebo
    res75 <- run_recession(cells, 1975)

    # Collect non-NULL models for etable
    models_list <- list(res83[[1]], res83[[2]], res75[[1]], res75[[2]])
    models_3way <- list(res83[[3]], res83[[4]], res75[[3]], res75[[4]])
    models_3way <- models_3way[!sapply(models_3way, is.null)]

    recession_dict <- c(mean_log_sls = "Mean log(Sales)",
                         treated_recession = "Treated $\\times$ Recession",
                         treated_single = "Treated $\\times$ Single-estab.",
                         single_recession = "Single-estab. $\\times$ Recession",
                         treated_recession_single = "Treated $\\times$ Recession $\\times$ Single-estab.",
                         exit_rate = "Exit rate",
                         fe = "State $\\times$ Industry $\\times$ Year FE")
    tex_out <- capture.output(etable(res83[[1]], res83[[2]], res75[[1]], res75[[2]],
                                     tex = TRUE, dict = recession_dict))
    writeLines(tex_out, file.path(OUTPUT_TAB, "recession.tex"))
    cat("    Created recession.tex\n")

    # Output 3-way FE version if available
    if (length(models_3way) > 0) {
      tex_3way <- capture.output(do.call(etable, c(models_3way, list(tex = TRUE, dict = recession_dict))))
      writeLines(tex_3way, file.path(OUTPUT_TAB, "recession_3way_fe.tex"))
      cat("    Created recession_3way_fe.tex\n")
    }

    rm(cells, res83, res75); gc()
  }, error = function(e) {
    cat("    Recession analysis failed:", e$message, "\n")
  })
}

cat("Marquette DiD analysis complete.\n")
