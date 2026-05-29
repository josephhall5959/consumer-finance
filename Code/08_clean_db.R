# 08_clean_db.R — Clean Dun & Bradstreet data
# retail.csv is the pre-processed analysis-ready retail establishment panel.
# db_70-85.csv contains raw D&B establishment records (partial: 1970 only).

cat("Cleaning D&B data...\n")

retail_file <- file.path(DATA_RAW, "DB", "retail.csv")
db_file <- file.path(DATA_RAW, "DB", "db_70-85.csv")

# --- retail.csv: pre-processed retail establishment panel ---
# retail.csv is large (5GB+). Rather than copying it, we create a symlink
# to avoid doubling disk usage. The high_acceptance column will be computed
# on-the-fly in 13_did_marquette.R.
if (file.exists(retail_file) && file.info(retail_file)$size > 1000) {
  cat("  Found retail.csv (", format(file.info(retail_file)$size, big.mark = ","),
      " bytes).\n")

  # Symlink to Generated/ so downstream scripts find it
  gen_retail <- file.path(DATA_GEN, "retail.csv")
  if (!file.exists(gen_retail)) {
    file.symlink(retail_file, gen_retail)
    cat("  Symlinked retail.csv to Data/Generated/\n")
  }

  # Create state-level aggregates (lightweight — reads only needed columns)
  cat("  Creating state-level aggregates from retail.csv...\n")
  retail_agg <- fread(retail_file,
                      select = c("YEAR", "STATE", "SLS", "entry", "exit", "card", "single"))
  state_agg <- retail_agg[, .(
    n = .N,
    sales = mean(as.numeric(SLS), na.rm = TRUE),
    entry = mean(entry, na.rm = TRUE),
    exit = mean(exit, na.rm = TRUE),
    n_card = sum(card, na.rm = TRUE),
    n_single = sum(single, na.rm = TRUE)
  ), by = .(YEAR, STATE)]
  fwrite(state_agg, file.path(DATA_GEN, "db_state.csv"))
  cat("  Wrote db_state.csv (state-level aggregates)\n")
  rm(retail_agg, state_agg)
  gc()

} else {
  cat("  retail.csv not found in Data/Raw/DB/\n")

  # Fallback: copy pre-generated state-level aggregates
  db_state_src <- file.path(DATA_RAW, "DB", "db_state.csv")
  if (file.exists(db_state_src)) {
    file.copy(db_state_src, file.path(DATA_GEN, "db_state.csv"), overwrite = TRUE)
    cat("  Copied db_state.csv (state-level aggregates)\n")
  }
}

# --- db_70-85.csv: raw D&B establishment records ---
if (file.exists(db_file) && file.info(db_file)$size > 1000) {
  cat("  Found db_70-85.csv (partial: 1970 only). Copying to Generated/\n")
  file.copy(db_file, file.path(DATA_GEN, "db_70-85.csv"), overwrite = TRUE)
} else {
  cat("  db_70-85.csv not available.\n")
}

# Copy other D&B intermediate files
for (f in c("uscities.csv", "fips_codes.csv", "db_state.csv")) {
  src <- file.path(DATA_RAW, "DB", f)
  if (file.exists(src)) file.copy(src, file.path(DATA_GEN, f), overwrite = TRUE)
}

cat("D&B cleaning complete.\n")
