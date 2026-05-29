# run_all.R — Master script for replication package
# Runs each script in a separate R process to keep memory bounded.
# Each subprocess inherits ROOT via environment variable.

cat("========================================\n")
cat("Credit Cards and Retail Firms\n")
cat("Replication Package\n")
cat("========================================\n\n")

# Determine Replication root
tryCatch({
  ROOT <- dirname(dirname(sys.frame(1)$ofile))
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) {
    ROOT <<- dirname(dirname(script_path))
  } else {
    ROOT <<- getwd()
  }
})
ROOT <- normalizePath(ROOT)

# ---------------------------------------------------------------------------
# Helper: run a script in a fresh R process
# ---------------------------------------------------------------------------
run_script <- function(script, description = NULL) {
  label <- if (!is.null(description)) description else script
  cat(sprintf("  [%s] %s ...\n", format(Sys.time(), "%H:%M:%S"), label))
  t0 <- proc.time()

  # Build the Rscript command.  Each child sources 00_setup.R (quietly), then the script.
  cmd <- sprintf(
    'Rscript -e "setwd(\\"%s\\"); suppressPackageStartupMessages(source(\\"Code/00_setup.R\\")); source(\\"%s\\")" 2>&1',
    gsub('"', '\\\\"', ROOT), gsub('"', '\\\\"', script)
  )
  rc <- system(cmd)

  elapsed <- round((proc.time() - t0)[["elapsed"]])
  if (rc != 0) {
    cat(sprintf("  *** FAILED: %s  (exit %d, %ds) ***\n", label, rc, elapsed))
    stop(sprintf("Script %s failed with exit code %d", script, rc))
  }
  cat(sprintf("  [%s] %s  (%ds)\n", format(Sys.time(), "%H:%M:%S"), label, elapsed))
  invisible(rc)
}

# ===== Phase 1: Data Cleaning ==============================================
cat("\n--- Phase 1: Data Cleaning ---\n")

run_script("Code/01_clean_usury.R",         "Usury laws")
run_script("Code/02_clean_scf.R",           "Survey of Consumer Finances")
run_script("Code/03_clean_cbp.R",           "County Business Patterns")
run_script("Code/04_clean_cps.R",           "Current Population Survey")
run_script("Code/05_clean_psid.R",          "PSID")
run_script("Code/06_clean_bea.R",           "BEA")
run_script("Code/07_clean_compustat.R",     "Compustat")
run_script("Code/08_clean_db.R",            "Dun & Bradstreet")
run_script("Code/09_pool_microdata.R",      "Pool microdata")

# ===== Phase 2: Analysis ===================================================
cat("\n--- Phase 2: Analysis ---\n")

run_script("Code/10_summary_stats.R",       "Summary statistics")
run_script("Code/11_did_bea_cbp.R",         "DiD: BEA/CBP")
run_script("Code/12_did_compustat.R",       "DiD: Compustat")
run_script("Code/13_did_marquette.R",       "DiD: Marquette")
run_script("Code/14_did_pooled.R",          "DiD: Pooled microdata")
# Code/15_analyze_scranton.R — superseded by 16_analyze_phonebooks.R
run_script("Code/16_analyze_phonebooks.R",  "Phonebook analysis")
run_script("Code/17_analyze_receivables.R", "Receivables")
run_script("Code/18_analyze_occupations.R", "Occupations")
run_script("Code/19_card_introduction.R",   "Card introduction dates")
run_script("Code/20_permutation_test.R",    "Permutation test")
run_script("Code/21_nested_logit.R",        "Nested logit model")
run_script("Code/22_macro_model.R",         "Macro model")

# ===== Done =================================================================
OUTPUT_FIG <- file.path(ROOT, "Output", "Figures")
OUTPUT_TAB <- file.path(ROOT, "Output", "Tables")

n_fig <- length(list.files(OUTPUT_FIG, pattern = "\\.(pdf|png)$"))
n_tab <- length(list.files(OUTPUT_TAB, pattern = "\\.tex$"))

cat("\n========================================\n")
cat("Pipeline complete!\n")
cat(sprintf("  Figures: %s  (%d files)\n", OUTPUT_FIG, n_fig))
cat(sprintf("  Tables:  %s  (%d files)\n", OUTPUT_TAB, n_tab))
cat("========================================\n")
