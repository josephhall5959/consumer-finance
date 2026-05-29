# 19_card_introduction.R — Card introduction timing analysis
# Source: "Card Introduction Summary Statistics .R"
# Uses Introduction Dates.xlsx from Yellow Pages data

cat("Analyzing card introduction...\n")

yp_dir <- file.path(DATA_RAW, "Yellow_Pages")

# --- 1. Phonebooks by year (from pages_include.csv) ---
pages_file <- file.path(yp_dir, "pages_include.csv")
if (file.exists(pages_file) && file.size(pages_file) > 0) {
  pages <- fread(pages_file)

  if ("year" %in% names(pages) || "Year" %in% names(pages)) {
    yr_col <- if ("year" %in% names(pages)) "year" else "Year"
    books_by_year <- pages[, .N, by = get(yr_col)]
    setnames(books_by_year, c("year", "n_books"))

    p <- ggplot(books_by_year, aes(x = year, y = n_books)) +
      geom_col(fill = "steelblue") +
      theme_classic() +
      labs(x = "Year", y = "Number of phonebooks")
    ggsave(file.path(OUTPUT_FIG, "books_by_year.pdf"), p, width = 9, height = 6)
    cat("  Created books_by_year.pdf\n")
  }
}

# --- 2. Card introduction dates (from Introduction Dates.xlsx) ---
intro_file <- file.path(yp_dir, "Introduction Dates.xlsx")
if (file.exists(intro_file) && file.size(intro_file) > 0) {
  cards.in <- tryCatch(read_excel(intro_file), error = function(e) {
    warning("Could not read Introduction Dates.xlsx: ", e$message)
    NULL
  })

  if (!is.null(cards.in) && "Card_intro" %in% names(cards.in)) {
    cards <- cards.in %>%
      mutate(
        precision = case_when(
          str_length(Card_intro) == 4 ~ "Exact",
          str_length(Card_intro) == 9 ~ "Range",
          str_length(Card_intro) == 5 & substr(Card_intro, 5, 5) == "+" ~ "Lower bound",
          str_length(Card_intro) == 5 & substr(Card_intro, 5, 5) == "-" ~ "Upper bound"
        ),
        card_year = as.numeric(substr(Card_intro, 1, 4))
      )

    # County-level lower bounds
    counties_lb <- cards %>%
      filter(!is.na(fipstate) & !is.na(fipscty)) %>%
      group_by(fipstate, fipscty) %>%
      summarise(intro_lb = min(card_year, na.rm = TRUE), .groups = "drop") %>%
      mutate(fips = paste0(str_pad(fipstate, 2, "left", "0"),
                           str_pad(fipscty, 3, "left", "0")))

    # Histogram of introduction lower bounds
    p <- ggplot(counties_lb, aes(x = intro_lb)) +
      geom_histogram(bins = 15, fill = "steelblue", alpha = 0.7, color = "white") +
      theme_classic() +
      labs(x = "Year of first credit card listing (lower bound)",
           y = "Number of counties")
    ggsave(file.path(OUTPUT_FIG, "intro_lb.png"), p, width = 9, height = 6, dpi = 300)
    cat("  Created intro_lb.png\n")

    # Map if sf is available
    tryCatch({
      shp_dir <- file.path(DATA_RAW, "Shapefiles")
      shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
      county_shp <- shp_files[grep("county|cb_.*_us_county", shp_files, ignore.case = TRUE)]

      if (length(county_shp) > 0) {
        counties_sf <- st_read(county_shp[1], quiet = TRUE)
        counties_sf$fips <- paste0(counties_sf$STATEFP, counties_sf$COUNTYFP)

        map_data <- left_join(counties_sf, counties_lb, by = "fips")

        p <- ggplot(map_data %>% filter(!is.na(intro_lb))) +
          geom_sf(aes(fill = intro_lb), linewidth = 0.05) +
          scale_fill_viridis_c(name = "Introduction year\n(lower bound)") +
          theme_void() +
          theme(legend.position = "bottom")
        ggsave(file.path(OUTPUT_FIG, "intro_lb_map.pdf"), p, width = 9, height = 6)
        cat("  Created intro_lb_map.pdf\n")
      }
    }, error = function(e) cat("  Map skipped:", e$message, "\n"))

  } else {
    cat("  Introduction Dates.xlsx missing Card_intro column. Skipping.\n")
  }
} else {
  cat("  Introduction Dates.xlsx not found. Skipping card introduction analysis.\n")
}

cat("Card introduction analysis complete.\n")
