# 00_setup.R — Package installation and path configuration

# Root path detection
if (requireNamespace("here", quietly = TRUE)) {
  ROOT <- here::here()
} else {
  # Fallback: assume working directory is Replication/
  ROOT <- getwd()
}

# Path definitions
DATA_RAW   <- file.path(ROOT, "Data", "Raw")
DATA_GEN   <- file.path(ROOT, "Data", "Generated")
OUTPUT_FIG <- file.path(ROOT, "Output", "Figures")
OUTPUT_TAB <- file.path(ROOT, "Output", "Tables")

# Create output directories if they don't exist
dir.create(DATA_GEN, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_TAB, showWarnings = FALSE, recursive = TRUE)

# Package installation helper
ensure_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    }
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}

# Required packages
PACKAGES <- c(
  "data.table", "dplyr", "tidyr", "stringr", "readr",
  "ggplot2", "gridExtra",
  "fixest", "modelsummary", "kableExtra",
  "haven", "readxl",
  "sf",
  "nleqslv", "numDeriv", "marginaleffects",
  "conflicted",
  "survey",
  "ipumsr"
)

# ggpubr is optional (provides theme_pubr); falls back to theme_classic
if (requireNamespace("ggpubr", quietly = TRUE)) {
  require(ggpubr)
} else {
  theme_pubr <- function(...) theme_classic(...)
}

ensure_packages(PACKAGES)

# Resolve common conflicts
suppressMessages({
  conflict_prefer("select", "dplyr", quiet = TRUE)
  conflict_prefer("filter", "dplyr", quiet = TRUE)
  conflict_prefer("lag", "dplyr", quiet = TRUE)
  conflict_prefer("between", "data.table", quiet = TRUE)
})

# Custom operator
'%!in%' <- function(x, y) !('%in%'(x, y))

cat("Setup complete. ROOT =", ROOT, "\n")
