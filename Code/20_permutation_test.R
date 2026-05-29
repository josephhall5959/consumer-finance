# 20_permutation_test.R — Permutation test for Compustat DiD

cat("Running permutation test...\n")

panel_file <- file.path(DATA_GEN, "cs_panel.csv")
reshaped_file <- file.path(DATA_GEN, "cs_db_reshaped.csv")

if (!file.exists(panel_file) || file.info(panel_file)$size < 100) {
  cat("  cs_panel.csv not found. Skipping permutation test.\n")
} else {
  panel <- fread(panel_file)

  usury <- fread(file.path(DATA_GEN, "usury.csv")) %>%
    dplyr::select(state, r_77) %>% distinct() %>%
    mutate(treatment_discrete = ifelse(r_77 < 18, 1, 0),
           treatment_continuous = (30 - r_77) / 100)
  usury$treatment_z <- (usury$treatment_continuous - mean(usury$treatment_continuous)) / sd(usury$treatment_continuous)

  # Load reshaped data if available
  if (file.exists(reshaped_file)) {
    reshaped <- fread(reshaped_file)
  } else {
    reshaped <- data.table(gvkey = integer(), sls_state = character(), sls_pct = numeric())
  }

  # Get unique states in the data
  states <- reshaped %>% dplyr::select(sls_state) %>%
    rbind(panel %>% dplyr::select(sls_state = state)) %>%
    mutate(sls_state = toupper(sls_state)) %>%
    distinct() %>%
    left_join(usury %>% dplyr::select(sls_state = state, treatment = treatment_discrete),
              by = "sls_state") %>%
    mutate(sls_state = tolower(sls_state)) %>%
    filter(!is.na(treatment))

  if (nrow(states) == 0 || !("rec_to_rev" %in% names(panel))) {
    cat("  Insufficient data for permutation test. Skipping.\n")
  } else {
    nperm <- 1000
    betas_cf <- rep(0, nperm)

    cat("  Running", nperm, "permutations...\n")

    for (i in 1:nperm) {
      if (i %% 100 == 0) cat("    Permutation", i, "\n")
      set.seed(i)
      vector_cf <- states %>% sample_frac(1, replace = TRUE) %>%
        dplyr::select(treatment_cf = treatment)

      treatment_cf <- states %>% mutate(treatment_cf = vector_cf$treatment_cf)

      treatment_cf_gvkey <- panel %>%
        dplyr::select(gvkey, state) %>% distinct() %>%
        left_join(reshaped, by = "gvkey", multiple = "all") %>%
        mutate(sls_state = ifelse(is.na(sls_state), tolower(state), sls_state),
               sls_pct = ifelse(is.na(sls_pct), 1, sls_pct)) %>%
        left_join(treatment_cf, by = "sls_state") %>%
        group_by(gvkey) %>%
        summarise(treatment_cf = weighted.mean(treatment_cf, sls_pct, na.rm = TRUE),
                  .groups = "drop")

      panel_cf <- panel %>% left_join(treatment_cf_gvkey, by = "gvkey")

      tryCatch({
        est_did <- feols(rec_to_rev ~ i(year, treatment_cf, 1978) | year + gvkey,
                         data = panel_cf %>% filter(year <= 1986))
        # Get the 1983 coefficient (5 years post)
        coefs <- summary(est_did)$coefficients
        target_coef <- coefs[grep("1983", names(coefs))]
        if (length(target_coef) > 0) {
          betas_cf[i] <- target_coef[1]
        }
      }, error = function(e) NULL)
    }

    # Plot permutation test results
    actual_coef <- -0.036  # From the main specification

    plot_data <- data.frame(betas_cf = betas_cf) %>%
      filter(betas_cf < 0.1 & betas_cf > -0.1)

    p <- ggplot(plot_data, aes(x = betas_cf)) +
      stat_ecdf(geom = "step") +
      geom_vline(xintercept = actual_coef, color = "red") +
      theme_classic() +
      labs(x = "Coefficient", y = "Empirical CDF")
    ggsave(file.path(OUTPUT_FIG, "permutation_test.pdf"), p, width = 9, height = 6)
    cat("  Created permutation_test.pdf\n")
  }
}

cat("Permutation test complete.\n")
