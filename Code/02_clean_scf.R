# 02_clean_scf.R — Load pre-cleaned Survey of Consumer Finances data
# Original cleaning was done in Stata ("01. clean SCF.do") which produces
# clean_ts.csv with state identifiers, usury rate merges, and treatment variables.
# The R translation cannot fully replicate this because it requires restricted-use
# geographic identifiers. We use the pre-cleaned output directly.

cat("02_clean_scf.R: Starting...\n")

scf_clean <- file.path(DATA_RAW, "SCF", "clean_ts.csv")

if (!file.exists(scf_clean)) {
  stop("Pre-cleaned SCF data not found at: ", scf_clean)
}

scf <- fread(scf_clean)
cat(sprintf("  Loaded %d observations across %d years\n",
            nrow(scf), length(unique(scf$year))))

# Verify key columns exist
required_cols <- c("weight", "income", "year", "has_card", "has_bank_card",
                   "has_store_card", "card_spend", "card_bal", "state", "r_77")
missing <- setdiff(required_cols, names(scf))
if (length(missing) > 0) {
  stop("clean_ts.csv is missing required columns: ", paste(missing, collapse = ", "))
}

# Copy to Generated
fwrite(scf, file.path(DATA_GEN, "clean_ts.csv"))
cat("  Copied to:", file.path(DATA_GEN, "clean_ts.csv"), "\n")

rm(scf); gc()
cat("02_clean_scf.R: Done.\n")
