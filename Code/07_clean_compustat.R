# 07_clean_compustat.R — Clean Compustat data

cat("Cleaning Compustat data...\n")

# Try to read fundamentals data
fund_file <- file.path(DATA_RAW, "Compustat", "fundamentals.csv")
public_file <- file.path(DATA_RAW, "Compustat", "public.csv")

if (file.exists(public_file) && file.info(public_file)$size > 100) {
  # Use pre-cleaned public.csv from Yens
  cat("  Using pre-cleaned public.csv\n")
  public_final <- fread(public_file)

  # Also copy state-level if available
  ps_file <- file.path(DATA_RAW, "Compustat", "public_state.csv")
  if (file.exists(ps_file)) {
    file.copy(ps_file, file.path(DATA_GEN, "public_state.csv"), overwrite = TRUE)
  }

  fwrite(public_final, file.path(DATA_GEN, "public.csv"))

} else if (file.exists(fund_file) && file.info(fund_file)$size > 100) {
  data.in <- fread(fund_file)

  # Create panel of firms which existed in 1978
  to_keep <- data.in %>%
    filter(!is.na(state)) %>%
    group_by(gvkey, state) %>%
    summarise(ymin = min(fyearq), ymax = max(fyearq), .groups = "drop") %>%
    filter(ymin < 1978 & ymax > 1978) %>%
    dplyr::select(gvkey, state)

  public_final <- to_keep %>%
    left_join(data.in, by = c("gvkey", "state")) %>%
    mutate(rec_to_assets = rectq / atq,
           rec_to_rev = rectq / revtq) %>%
    filter(revtq > 0 & atq > 0 & rectq >= 0 & rec_to_assets < 1 & rec_to_rev < 1)

  fwrite(public_final, file.path(DATA_GEN, "public.csv"))

  # State-level summaries
  state_all <- public_final %>%
    group_by(state, datadate) %>%
    summarise(revenue = mean(revtq), assets = mean(atq),
              rec_to_assets = mean(rec_to_assets), rec_to_rev = mean(rec_to_rev),
              n_public = n(), .groups = "drop")

  state_retail <- public_final %>%
    filter(sic >= 5000 & sic < 6000) %>%
    group_by(state, datadate) %>%
    summarise(revenue_retail = mean(revtq), assets_retail = mean(atq),
              rec_to_assets_retail = mean(rec_to_assets), rec_to_rev_retail = mean(rec_to_rev),
              n_public_retail = n(), .groups = "drop")

  public_state <- state_all %>%
    left_join(state_retail, by = c("state", "datadate"))

  fwrite(public_state, file.path(DATA_GEN, "public_state.csv"))
} else {
  warning("No Compustat fundamentals data found. Using available intermediate files.")
  # Copy any available intermediate files
  for (f in c("cs_panel.csv", "cs_xsection.csv", "cs_panel_qly.csv", "cs_profit_growth.csv")) {
    src <- file.path(DATA_RAW, "Compustat", f)
    if (file.exists(src)) file.copy(src, file.path(DATA_GEN, f), overwrite = TRUE)
  }
}

# Copy other intermediate Compustat files that downstream scripts need
for (f in c("cs_panel.csv", "cs_retail_all_matched.csv", "cs_db_reshaped.csv", "cs_with_db_data_v1.csv", "retail_gvkeys.csv")) {
  src <- file.path(DATA_RAW, "Compustat", f)
  if (file.exists(src) && !file.exists(file.path(DATA_GEN, f))) {
    file.copy(src, file.path(DATA_GEN, f), overwrite = TRUE)
  }
}

cat("Compustat cleaning complete.\n")
