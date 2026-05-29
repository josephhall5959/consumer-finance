# 01_clean_usury.R — Clean usury rate data and create related outputs
# Translated from Stata: "00. clean.do"
#
# Assumes source("Code/00_setup.R") has been run, providing:
#   ROOT, DATA_RAW, DATA_GEN, OUTPUT_FIG, OUTPUT_TAB
#   and all required packages loaded (data.table, dplyr, tidyr, stringr,
#   readr, ggplot2, fixest, haven, readxl, sf, kableExtra, conflicted)

cat("01_clean_usury.R: Starting...\n")

# =============================================================================
# 1. Cross-country credit card ownership regressions
# =============================================================================

cat("  Reading cross_country.xlsx...\n")
cc <- readxl::read_excel(file.path(DATA_RAW, "cross_country.xlsx"),
                         sheet = "Sheet1")

# Regression: CreditCardOwnership ~ Year + Country FE, clustered by Country
# (xi:reg with i.Country in Stata => fixed effects for Country)
if ("CreditCardOwnership" %in% names(cc) && "Year" %in% names(cc) &&
    "Country" %in% names(cc)) {

  cc_reg <- cc[!is.na(cc$CreditCardOwnership), ]
  if (nrow(cc_reg) > 0) {
    reg1 <- fixest::feols(CreditCardOwnership ~ Year | Country,
                          data = cc_reg, cluster = ~Country)
    cat("  Regression 1 (CreditCardOwnership ~ Year | Country):\n")
    print(summary(reg1))
  }

  # Regression: StoreCreditBorrowing ~ Year + Country FE, clustered by Country
  cc_reg2 <- cc[!is.na(cc$StoreCreditBorrowing), ]
  if ("StoreCreditBorrowing" %in% names(cc) && nrow(cc_reg2) > 0) {
    reg2 <- fixest::feols(StoreCreditBorrowing ~ Year | Country,
                          data = cc_reg2, cluster = ~Country)
    cat("  Regression 2 (StoreCreditBorrowing ~ Year | Country):\n")
    print(summary(reg2))
  }

  # Keep only obs with both variables non-missing, then keep latest year per country
  cc_both <- cc %>%
    dplyr::filter(!is.na(CreditCardOwnership) & !is.na(StoreCreditBorrowing)) %>%
    dplyr::arrange(Country, desc(Year)) %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::filter(rank == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-rank)
}

# =============================================================================
# 2. Read and reshape usury rates: wide -> long, handle "nolimit"
# =============================================================================

cat("  Reading usury_rates.csv...\n")
usury_wide <- readr::read_csv(file.path(DATA_RAW, "usury_rates.csv"),
                               show_col_types = FALSE)

# Rename first column to 'state' to match Stata variable naming
names(usury_wide)[names(usury_wide) == "State"] <- "state"

# Reshape from wide to long: columns y1971..y1990 -> (year, r_usury)
usury <- usury_wide %>%
  tidyr::pivot_longer(
    cols = starts_with("y"),
    names_to = "year",
    names_prefix = "y",
    values_to = "r_usury"
  ) %>%
  dplyr::mutate(year = as.integer(year))

# Create nolimit indicator (where r_usury == "nolimit")
usury <- usury %>%
  dplyr::mutate(
    nolimit = as.integer(r_usury == "nolimit"),
    # Force convert to numeric (non-numeric strings like "nolimit" become NA)
    r_usury = suppressWarnings(as.numeric(r_usury))
  )

# For nolimit states: set rate to max(r_usury) + 5 within each year
usury <- usury %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(
    rmax = max(r_usury, na.rm = TRUE),
    r_usury = ifelse(nolimit == 1, rmax + 5, r_usury)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-rmax)

# Drop observations with missing r_usury
usury <- usury %>%
  dplyr::filter(!is.na(r_usury))

cat("  Usury data: ", nrow(usury), " observations\n")

# =============================================================================
# 3. Merge with states.csv
# =============================================================================

cat("  Reading states.csv and merging...\n")
states <- readr::read_csv(file.path(DATA_RAW, "states.csv"),
                           show_col_types = FALSE)

# Stata code: ren state NAME; drop abbrev; ren code state
# states.csv columns: State, Abbrev, Code, st, state_psid, msv_index
# We need: NAME (= State full name), state (= Code = 2-letter abbreviation)
states_clean <- states %>%
  dplyr::rename(NAME = State) %>%
  dplyr::select(-Abbrev) %>%
  dplyr::rename(state_code = Code)

# Merge: states (1) to usury (n) by state code
# In usury, 'state' is the 2-letter code; in states_clean, 'state_code' is the 2-letter code
usury <- usury %>%
  dplyr::inner_join(states_clean, by = c("state" = "state_code"))

cat("  After merge with states: ", nrow(usury), " observations\n")

# =============================================================================
# 4. Read Fed Funds rate and inflation data
# =============================================================================

cat("  Reading FEDFUNDS.csv...\n")
fedfunds <- readr::read_csv(file.path(DATA_RAW, "FEDFUNDS.csv"),
                              show_col_types = FALSE)
fedfunds <- fedfunds %>%
  dplyr::mutate(year = as.integer(substr(DATE, 1, 4))) %>%
  dplyr::rename(fedfunds = FEDFUNDS) %>%
  dplyr::select(year, fedfunds)

cat("  Reading FPCPITOTLZGUSA.csv...\n")
inflation <- readr::read_csv(file.path(DATA_RAW, "FPCPITOTLZGUSA.csv"),
                               show_col_types = FALSE)
inflation <- inflation %>%
  dplyr::mutate(year = as.integer(substr(DATE, 1, 4))) %>%
  dplyr::rename(inflation = FPCPITOTLZGUSA) %>%
  dplyr::select(year, inflation)

# =============================================================================
# 5. Create time series: average usury rate + nolimit count by year,
#    merged with fed funds and inflation
# =============================================================================

cat("  Creating time series...\n")
usury_ts <- usury %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    r_usury = mean(r_usury, na.rm = TRUE),
    nolimit = sum(nolimit, na.rm = TRUE),
    .groups = "drop"
  )

usury_ts <- usury_ts %>%
  dplyr::left_join(fedfunds, by = "year") %>%
  dplyr::left_join(inflation, by = "year")

# =============================================================================
# 6. Create 1977 cross-section and 1990 cross-section, then combine
# =============================================================================

cat("  Creating usury cross-section (1977 & 1990)...\n")

# 1977 cross-section
u_77 <- usury %>%
  dplyr::filter(year == 1977) %>%
  dplyr::arrange(state) %>%
  dplyr::rename(r_77 = r_usury) %>%
  dplyr::select(NAME, state, r_77)

# 1990 cross-section
u_90 <- usury %>%
  dplyr::filter(year == 1990) %>%
  dplyr::rename(r_90 = r_usury) %>%
  dplyr::select(NAME, state, r_90)

# Merge 1977 and 1990 to create cross-section with change variable
usury_xsection <- u_90 %>%
  dplyr::inner_join(u_77, by = c("NAME", "state")) %>%
  dplyr::mutate(r_change = r_90 - r_77)

# Merge cross-section variables back to usury panel
usury <- usury %>%
  dplyr::left_join(
    usury_xsection %>% dplyr::select(state, r_77, r_90, r_change),
    by = "state"
  )

# =============================================================================
# 7. Save outputs to DATA_GEN
# =============================================================================

cat("  Saving usury.csv...\n")
readr::write_csv(usury, file.path(DATA_GEN, "usury.csv"))

cat("  Saving usury_xsection.csv...\n")
readr::write_csv(usury_xsection, file.path(DATA_GEN, "usury_xsection.csv"))

# =============================================================================
# 8. Create choropleth map of 1977 usury rates using sf + ggplot2
# =============================================================================

cat("  Creating choropleth map (rate_77.png)...\n")

tryCatch({
  # Try to read shapefiles from DATA_RAW/Shapefiles/
  shapefile_path <- file.path(DATA_RAW, "Shapefiles", "Shapefiles", "st99_d90.shp")

  if (!file.exists(shapefile_path) || file.size(shapefile_path) == 0) {
    # Try alternative: use US Census tigris or a bundled shapefile
    # Fall back to spData or maps package if available
    if (requireNamespace("maps", quietly = TRUE)) {
      cat("    Shapefiles not available; using maps package fallback...\n")
      us_states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
      us_states <- us_states %>%
        dplyr::mutate(NAME = stringr::str_to_title(ID))
    } else {
      stop("No shapefiles or maps package available for choropleth.")
    }
  } else {
    us_states <- sf::st_read(shapefile_path, quiet = TRUE)
  }

  # Merge map data with 1977 cross-section
  # Exclude Alaska and Hawaii (as in Stata code)
  map_data <- us_states %>%
    dplyr::filter(!NAME %in% c("Alaska", "Hawaii"))

  map_data <- map_data %>%
    dplyr::left_join(
      usury_xsection %>% dplyr::select(NAME, r_77),
      by = "NAME"
    )

  # Fill missing r_77 with 18 (as in Stata: replace r_77 = 18 if missing(r_77))
  map_data <- map_data %>%
    dplyr::mutate(r_77 = ifelse(is.na(r_77), 18, r_77))

  # Create factor for discrete color mapping (clm(unique) in Stata's spmap)
  map_data <- map_data %>%
    dplyr::mutate(r_77_factor = factor(r_77))

  # Plot choropleth
  p <- ggplot2::ggplot(map_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = r_77_factor), color = "white", size = 0.2) +
    ggplot2::scale_fill_brewer(palette = "Blues", name = "Usury Rate (1977)") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::ggtitle("Usury Rate Limits, 1977")

  ggplot2::ggsave(file.path(OUTPUT_FIG, "rate_77.png"), plot = p,
                  width = 10, height = 6, dpi = 300)
  cat("    Saved rate_77.png\n")

}, error = function(e) {
  cat("    WARNING: Could not generate choropleth map.\n")
  cat("    Reason:", conditionMessage(e), "\n")
  cat("    Skipping map generation.\n")
})

cat("01_clean_usury.R: Done.\n")
