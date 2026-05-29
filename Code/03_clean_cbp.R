# 03_clean_cbp.R — Clean County Business Patterns data
# Translated from Stata: "03. clean CBP.do"
#
# Assumes source("Code/00_setup.R") has been run, providing:
#   ROOT, DATA_RAW, DATA_GEN, OUTPUT_FIG, OUTPUT_TAB
#   and all required packages loaded (data.table, dplyr, tidyr, stringr,
#   readr, ggplot2, fixest, haven, readxl, sf, kableExtra, conflicted)

cat("03_clean_cbp.R: Starting...\n")

# Base path for CBP data
CBP_DIR <- file.path(DATA_RAW, "CBP")

# =============================================================================
# 1. Check what data files are available
# =============================================================================

# CBP data comes in fixed-width text files (cbpYYst.txt) for state-level data.
# The format varies by decade:
#   - 1970s-1980s: older fixed-width format with SIC codes
#   - 1990s-1997: state-level with SIC codes, slightly different layout
#   - 1998+: NAICS codes (different industry classification)
#
# If pre-processed CSV files are available, prefer those.

cat("  Checking available CBP data...\n")

# Check for pre-processed files first
sic2_csv_path <- file.path(CBP_DIR, "sic2.csv")
panel_csv_path <- file.path(CBP_DIR, "panel.csv")
retail_csv_path <- file.path(CBP_DIR, "retail.csv")

has_sic2_csv   <- file.exists(sic2_csv_path) && file.size(sic2_csv_path) > 0
has_panel_csv  <- file.exists(panel_csv_path) && file.size(panel_csv_path) > 0
has_retail_csv <- file.exists(retail_csv_path) && file.size(retail_csv_path) > 0

# =============================================================================
# 2. If pre-processed CSV is available, read it directly
# =============================================================================

if (has_sic2_csv) {
  cat("  Found pre-processed sic2.csv — reading directly.\n")
  cbp_sic2 <- readr::read_csv(sic2_csv_path, show_col_types = FALSE)
  cat(sprintf("  Loaded %d rows from sic2.csv\n", nrow(cbp_sic2)))

} else if (has_panel_csv) {
  cat("  Found pre-processed panel.csv — reading and filtering to retail SIC.\n")
  cbp_panel <- readr::read_csv(panel_csv_path, show_col_types = FALSE)

  # Filter to retail SIC codes (52-59) and aggregate to 2-digit level
  sic_col <- intersect(c("sic", "SIC", "sic2", "industry"), names(cbp_panel))
  if (length(sic_col) > 0) {
    cbp_sic2 <- cbp_panel %>%
      dplyr::filter(.data[[sic_col[1]]] >= 52 & .data[[sic_col[1]]] <= 59)
  } else {
    cbp_sic2 <- cbp_panel
  }

} else {
  # =============================================================================
  # 3. Read raw fixed-width text files year by year
  # =============================================================================

  cat("  No pre-processed CSV found. Attempting to read raw CBP text files...\n")

  # Determine which years have state-level data files
  cbp_years <- list.dirs(CBP_DIR, full.names = FALSE, recursive = FALSE)
  cbp_years <- cbp_years[grepl("^\\d{4}$", cbp_years)]
  cbp_years <- sort(as.integer(cbp_years))

  cat(sprintf("  Found year directories: %s\n", paste(cbp_years, collapse = ", ")))

  # -------------------------------------------------------------------------
  # Define reading functions for different CBP file formats
  # -------------------------------------------------------------------------

  # CBP state-level files (cbpYYst.txt) use fixed-width format.
  # The column layout changed over time. The key columns are:
  #   fipstate (2 chars) — state FIPS code
  #   sic (4 chars)      — SIC industry code
  #   emp / est          — employment and establishment counts
  #   emp_nf             — employment noise flag
  #   Various employment-size buckets (n1_4, n5_9, n10_19, etc.)

  read_cbp_state <- function(year) {
    # Construct the file path
    yy <- substr(as.character(year), 3, 4)
    filepath <- file.path(CBP_DIR, as.character(year), paste0("cbp", yy, "st.txt"))

    if (!file.exists(filepath)) {
      warning(sprintf("  CBP state file not found for %d: %s", year, filepath))
      return(NULL)
    }

    if (file.size(filepath) == 0) {
      warning(sprintf("  CBP state file is empty (stub) for %d: %s", year, filepath))
      return(NULL)
    }

    cat(sprintf("    Reading %d: %s\n", year, basename(filepath)))

    # Try reading as CSV first (some CBP files are comma-delimited despite .txt)
    df <- tryCatch({
      readr::read_csv(filepath, show_col_types = FALSE)
    }, error = function(e) {
      NULL
    })

    # If CSV read worked and has reasonable columns, use it
    if (!is.null(df) && ncol(df) > 3) {
      names(df) <- tolower(names(df))
      df$year <- year
      return(df)
    }

    # Try fixed-width format
    # The CBP fixed-width layout varies by decade. Common layouts:
    #
    # Pre-1986 format (approx):
    #   Col 1-2:   FIPS state code
    #   Col 3-6:   SIC code
    #   Col 7-15:  Total establishments
    #   Col 16-24: Total mid-March employment
    #   Col 25-33: First quarter payroll ($1000)
    #   Col 34-42: Annual payroll ($1000)
    #   ... employment size class columns follow
    #
    # Post-1986 format is similar but may have additional columns.

    df <- tryCatch({
      # Attempt a generic fixed-width read with common column positions
      # These positions are approximate and work for many years
      if (year <= 1985) {
        fwf_positions <- readr::fwf_widths(
          widths = c(2, 4, 1, 9, 9, 12, 12,
                     9, 9, 9, 9, 9, 9, 9, 9, 9, 9),
          col_names = c("fipstate", "sic", "emp_flag",
                        "est", "emp", "qp1", "ap",
                        "n1_4", "n5_9", "n10_19", "n20_49",
                        "n50_99", "n100_249", "n250_499",
                        "n500_999", "n1000_plus", "n1000_1")
        )
      } else {
        fwf_positions <- readr::fwf_widths(
          widths = c(2, 4, 1, 9, 1, 9, 12, 1, 12,
                     9, 9, 9, 9, 9, 9, 9, 9, 9, 9),
          col_names = c("fipstate", "sic", "emp_flag",
                        "est", "emp_nf", "emp", "qp1", "qp1_nf", "ap",
                        "n1_4", "n5_9", "n10_19", "n20_49",
                        "n50_99", "n100_249", "n250_499",
                        "n500_999", "n1000_plus", "n1000_1")
        )
      }

      readr::read_fwf(filepath, col_positions = fwf_positions,
                      show_col_types = FALSE)
    }, error = function(e) {
      warning(sprintf("    Could not parse fixed-width file for %d: %s",
                      year, conditionMessage(e)))
      NULL
    })

    if (!is.null(df)) {
      names(df) <- tolower(names(df))
      df$year <- year
    }

    return(df)
  }

  # -------------------------------------------------------------------------
  # Read all available years
  # -------------------------------------------------------------------------

  all_cbp <- list()

  for (yr in cbp_years) {
    result <- tryCatch(
      read_cbp_state(yr),
      error = function(e) {
        warning(sprintf("  Error reading CBP %d: %s", yr, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(result)) {
      all_cbp[[as.character(yr)]] <- result
    }
  }

  if (length(all_cbp) == 0) {
    cat("  WARNING: No CBP data files could be read.\n")
    cat("  This may be because the raw data files are placeholder stubs (0 bytes).\n")
    cat("  To run this script, place the actual CBP text files in:\n")
    cat("    ", CBP_DIR, "\n")
    cat("  Or provide a pre-processed sic2.csv, panel.csv, or retail.csv.\n")
  } else {
    # -----------------------------------------------------------------------
    # Combine all years
    # -----------------------------------------------------------------------

    cat("  Combining CBP data across years...\n")
    cbp_all <- dplyr::bind_rows(all_cbp)

    # Standardize column names
    names(cbp_all) <- tolower(names(cbp_all))

    # -----------------------------------------------------------------------
    # 4. Filter to retail SIC codes (52-59)
    # -----------------------------------------------------------------------

    cat("  Filtering to retail SIC codes (52-59)...\n")

    # SIC codes may be character or numeric; ensure consistent handling
    if ("sic" %in% names(cbp_all)) {
      cbp_all <- cbp_all %>%
        dplyr::mutate(
          sic_clean = stringr::str_trim(as.character(sic)),
          # Extract 2-digit SIC code
          sic2 = as.integer(substr(sic_clean, 1, 2))
        )

      # Keep only retail SIC codes: 52 (Building materials), 53 (General merch),
      # 54 (Food stores), 55 (Auto dealers), 56 (Apparel), 57 (Furniture),
      # 58 (Eating/drinking), 59 (Misc retail)
      cbp_retail <- cbp_all %>%
        dplyr::filter(sic2 >= 52 & sic2 <= 59)
    } else {
      # If no SIC column, try industry/naics columns
      cat("  WARNING: No 'sic' column found. Keeping all data.\n")
      cbp_retail <- cbp_all
      cbp_retail$sic2 <- NA_integer_
    }

    # -----------------------------------------------------------------------
    # 5. Aggregate to 2-digit SIC level by state and year
    # -----------------------------------------------------------------------

    cat("  Aggregating to 2-digit SIC by state and year...\n")

    # Identify the state column
    state_col <- intersect(c("fipstate", "fips_state", "state", "st"),
                           names(cbp_retail))

    if (length(state_col) == 0) {
      cat("  WARNING: No state identifier column found. Using all data as-is.\n")
      cbp_retail$fipstate <- NA
      state_col <- "fipstate"
    }

    # Ensure numeric columns are numeric
    num_cols <- c("est", "emp", "ap", "qp1",
                  "n1_4", "n5_9", "n10_19", "n20_49",
                  "n50_99", "n100_249", "n250_499",
                  "n500_999", "n1000_plus")
    for (col in intersect(num_cols, names(cbp_retail))) {
      cbp_retail[[col]] <- suppressWarnings(as.numeric(cbp_retail[[col]]))
    }

    # Aggregate: sum establishments and employment by state x year x SIC-2
    agg_vars <- intersect(num_cols, names(cbp_retail))

    if (length(agg_vars) > 0) {
      cbp_sic2 <- cbp_retail %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(state_col[1], "year", "sic2")))) %>%
        dplyr::summarise(
          dplyr::across(dplyr::all_of(agg_vars), ~ sum(.x, na.rm = TRUE)),
          .groups = "drop"
        )
    } else {
      cbp_sic2 <- cbp_retail %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(state_col[1], "year", "sic2")))) %>%
        dplyr::summarise(n_obs = dplyr::n(), .groups = "drop")
    }
  }
}

# =============================================================================
# 6. Save output to DATA_GEN
# =============================================================================

if (exists("cbp_sic2") && nrow(cbp_sic2) > 0) {
  cat(sprintf("  Final dataset: %d rows\n", nrow(cbp_sic2)))
  cat("  Saving cbp_sic2.csv...\n")
  readr::write_csv(cbp_sic2, file.path(DATA_GEN, "cbp_sic2.csv"))
  cat("  Saved to:", file.path(DATA_GEN, "cbp_sic2.csv"), "\n")
} else {
  cat("  WARNING: No CBP data was produced.\n")
  cat("  Ensure that either:\n")
  cat("    (a) Raw CBP text files are in ", CBP_DIR, "/<year>/cbpYYst.txt\n")
  cat("    (b) Pre-processed CSV files (sic2.csv or panel.csv) are in ", CBP_DIR, "\n")
}

cat("03_clean_cbp.R: Done.\n")
