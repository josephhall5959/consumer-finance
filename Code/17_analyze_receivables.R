# 17_analyze_receivables.R — Analyze firm receivables (QFR data)

cat("Analyzing receivables...\n")

qfr_file <- file.path(DATA_RAW, "QFR", "qfr_all.csv")
if (!file.exists(qfr_file)) {
  stop("QFR data not found at: ", qfr_file)
}

qfr <- fread(qfr_file)
cat(sprintf("  Loaded %d rows from qfr_all.csv\n", nrow(qfr)))

# Extract numeric industry code
qfr[, industry_code := as.numeric(stringr::str_extract(`INDUSTRY CODE`, "\\d+"))]

# Classify into wholesale and retail trade
# SIC: Wholesale 50-51, Retail 53-54; NAICS: Wholesale 421-422, Retail 445-450
trade <- qfr %>%
  mutate(Industry = case_when(
    industry_code %in% c(50, 51, 421, 422) ~ "Wholesale trade",
    industry_code %in% c(53, 54, 445, 448, 450) ~ "Retail trade"
  )) %>%
  filter(!is.na(Industry) & TOTRECEV > 0) %>%
  group_by(Industry, year, quarter) %>%
  summarise(
    SALES = sum(as.numeric(SALES), na.rm = TRUE),
    TOTRECEV = sum(as.numeric(TOTRECEV), na.rm = TRUE),
    TOTASSET = sum(as.numeric(TOTASSET), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(rec_rev = TOTRECEV / SALES,
         date = year + 0.25 * quarter - 0.25)

if (nrow(trade) > 1) {
  p <- ggplot(trade %>% filter(year > 1981),
              aes(x = date, y = rec_rev, color = Industry, linetype = Industry)) +
    geom_line() +
    theme_classic() +
    scale_y_continuous(name = "Receivables/Revenue",
                       labels = scales::percent, limits = c(0, .5)) +
    scale_x_continuous(name = "Date (Quarterly)", limits = c(1980, 2020)) +
    labs(caption = "Source: Quarterly Financial Reports. Includes firms with >50m in assets")
  ggsave(file.path(OUTPUT_FIG, "receivables.pdf"), p, width = 9, height = 6)
  cat("  Created receivables.pdf\n")
}

rm(qfr, trade); gc()
cat("Receivables analysis complete.\n")
