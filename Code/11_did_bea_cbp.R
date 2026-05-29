# 11_did_bea_cbp.R — Difference-in-differences using BEA and CBP data
# Matches original "04. dif-in-dif.R"

cat("Running BEA/CBP DiD analysis...\n")

pretty_did <- function(did, ref_yr, ylab = "Outcome", title = "", breaks = 1) {
  name <- names(summary(did)$coefficients)
  beta <- summary(did)$coefficients
  se <- summary(did)$se
  toplot <- data.frame(name, beta, se) %>%
    mutate(year = readr::parse_number(name))
  zero <- data.frame(name = "_", year = ref_yr, beta = 0, se = 0)
  toplot <- rbind(toplot, zero)
  start <- min(toplot$year)
  end <- max(toplot$year)

  ggplot(data = toplot, aes(x = year, y = beta)) +
    geom_line() +
    geom_ribbon(aes(ymin = beta - 2 * se, ymax = beta + 2 * se), alpha = 0.2) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = ref_yr, linetype = "dashed") +
    scale_x_continuous(name = "Year", breaks = seq(start, end, breaks)) +
    scale_y_continuous(name = ylab) +
    ggtitle(title)
}

pretty_3did <- function(did, ref_yr, ylab = "Outcome", title = "", breaks = 1) {
  name <- names(summary(did)$coefficients)
  beta <- summary(did)$coefficients
  se <- summary(did)$se
  # First half of coefficients are the triple-diff (treated_clothing interaction)
  toplot <- data.frame(name, beta, se)[1:(length(name) / 2), ] %>%
    mutate(year = readr::parse_number(name))
  zero <- data.frame(name = "_", year = ref_yr, beta = 0, se = 0)
  toplot <- rbind(toplot, zero)
  start <- min(toplot$year)
  end <- max(toplot$year)

  ggplot(data = toplot, aes(x = year, y = beta)) +
    geom_line() +
    geom_ribbon(aes(ymin = beta - 2 * se, ymax = beta + 2 * se), alpha = 0.2) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = ref_yr, linetype = "dashed") +
    scale_x_continuous(name = "Year", breaks = seq(start, end, breaks)) +
    scale_y_continuous(name = ylab) +
    ggtitle(title)
}


# ============================================================
# PART 1: BEA State-Industry DiD (log levels)
# ============================================================

bea_file <- file.path(DATA_GEN, "bea.csv")
if (file.exists(bea_file) && file.info(bea_file)$size > 100) {
  cat("  Loading BEA data...\n")
  data <- fread(bea_file)

  data <- data[, .(GDP = sum(GDP, na.rm = TRUE), Employment = sum(Employment, na.rm = TRUE),
                    Wages = sum(Wages, na.rm = TRUE),
                    cpi_delneg = mean(cpi_delneg)),
               by = c("year", "fipstate", "tradable", "retail", "r_77", "treated")] %>%
    as_tibble() %>%
    group_by(fipstate, tradable, retail, r_77, treated) %>%
    arrange(fipstate, year) %>%
    mutate(GDP = GDP / cpi_delneg,
           Wages = Wages / cpi_delneg,
           wage_rate = Wages / Employment,
           gdp_growth = GDP / lag(GDP) - 1,
           emp_growth = Employment / lag(Employment) - 1,
           wage_growth = wage_rate / lag(wage_rate) - 1,
           price_growth = cpi_delneg / lag(cpi_delneg) - 1) %>%
    as.data.table()

  data[, treated_tradable := treated * tradable]
  data[, treated_nontradable := treated * (1 - tradable)]
  data[, treated_retail := treated * retail]
  data[, log_gdp := log(GDP)]
  data[, log_emp := log(Employment)]
  data[, log_wage := log(wage_rate)]
  data[, log_cpi := log(cpi_delneg)]

  # GDP figures
  est_did <- feols(log_gdp ~ i(year, treated, 1977) | year + fipstate,
                   cluster = ~fipstate, data = data[retail == 1])
  p <- pretty_did(est_did, 1977, "GDP (log points)", "Retail", breaks = 3)
  ggsave(file.path(OUTPUT_FIG, "gdp_retail.pdf"), p, width = 9, height = 6)

  est_did <- feols(log_gdp ~ i(year, treated_retail, 1977) + i(year, treated, 1977) | year + fipstate + retail + tradable,
                   cluster = ~fipstate, data = data)
  p <- pretty_3did(est_did, 1977, "GDP (log points)", "Retail vs. Other", breaks = 3)
  ggsave(file.path(OUTPUT_FIG, "gdp_retail_3d.pdf"), p, width = 9, height = 6)

  # Employment figures
  est_did <- feols(log_emp ~ i(year, treated, 1977) | year + fipstate,
                   cluster = ~fipstate, data = data[retail == 1])
  p <- pretty_did(est_did, 1977, "Employment (log points)", "Retail", breaks = 3)
  ggsave(file.path(OUTPUT_FIG, "emp_retail.pdf"), p, width = 9, height = 6)

  est_did <- feols(log_emp ~ i(year, treated_retail, 1977) + i(year, treated, 1977) | year + retail + tradable,
                   cluster = ~fipstate, data = data)
  p <- pretty_3did(est_did, 1977, "Employment (log points)", "Retail vs. Other", breaks = 3)
  ggsave(file.path(OUTPUT_FIG, "emp_retail_3d.pdf"), p, width = 9, height = 6)

  # Wage figures
  est_did <- feols(log_wage ~ i(year, treated, 1977) | year + fipstate,
                   cluster = ~fipstate, data = data[retail == 1])
  p <- pretty_did(est_did, 1977, "Wages per worker (log points)", "Retail", breaks = 3)
  ggsave(file.path(OUTPUT_FIG, "wage_retail.pdf"), p, width = 9, height = 6)

  est_did <- feols(log_wage ~ i(year, treated_retail, 1977) + i(year, treated, 1977) | year + retail + tradable,
                   cluster = ~fipstate, data = data)
  p <- pretty_3did(est_did, 1977, "Wages per worker (log points)", "Retail vs. Other", breaks = 3)
  ggsave(file.path(OUTPUT_FIG, "wage_retail_3d.pdf"), p, width = 9, height = 6)

  cat("  BEA figures complete.\n")
  rm(data); gc()
} else {
  cat("  BEA data not found, skipping.\n")
}


# ============================================================
# PART 2: CBP County-Level Triple-Difference (Clothing × Treated × Year)
# Matches original "04. dif-in-dif.R" lines 210-286
# ============================================================

cbp_file <- file.path(DATA_GEN, "cbp_sic2.csv")
sic2_map_file <- file.path(DATA_RAW, "CBP", "sic2_map.xlsx")

if (file.exists(cbp_file) && file.info(cbp_file)$size > 100) {
  cat("  Loading CBP county-level data...\n")
  data <- fread(cbp_file, select = c("stcode", "cocode", "y1977", "year", "sic2",
                                      "est", "n1_19", "n20_49", "n50_99",
                                      "n100_249", "n250_499", "n500"))
  # Filter to 1974-1989 (pre-1974 definitions change)
  data <- data[year >= 1974 & year < 1990]
  cat(sprintf("  %d rows after year filter\n", nrow(data)))

  # Unit = industry-state-county
  data[, unit := paste0(sic2, "_", stcode, "_", cocode)]

  # Treatment: states with usury rate of 12 or 15 in 1977
  data[, treated := ifelse(y1977 %in% c(12, 15), 1, 0)]

  # Recompute total establishments from size bins (matching original)
  data[, est := n1_19 + n20_49 + n50_99 + n100_249 + n250_499 + n500]
  data[, frac_small := n1_19 / est]

  # Merge sic2_map for tradable/retail/clothing flags
  if (file.exists(sic2_map_file)) {
    sic2_map <- readxl::read_excel(sic2_map_file) %>%
      dplyr::select(-N, -prop)
    data <- data %>% left_join(sic2_map, by = "sic2")
  } else {
    # Fallback: manually define retail and clothing
    cat("  sic2_map.xlsx not found, using manual SIC classification\n")
    data[, retail := ifelse(sic2 >= 52 & sic2 <= 59, 1, 0)]
    data[, clothing := ifelse(sic2 == 56, 1, 0)]
    data[, tradable := 0]  # approximate
  }

  # Interaction terms
  data[, treated_clothing := treated * clothing]
  data[, treated_tradable := treated * tradable]
  data[, treated_nontradable := treated * (1 - tradable)]

  cat("  Running triple-diff regressions...\n")

  # === Triple-Difference: Clothing × Treated × Year ===
  # est.pdf — Total establishments
  est_did <- feols(est ~ i(year, treated_clothing, 1977) + i(year, treated, 1977) | unit + year,
                   cluster = ~stcode,
                   data = data[year >= 1974 & year < 1984])
  p <- pretty_3did(est_did, 1977, "Count", "Establishments")
  ggsave(file.path(OUTPUT_FIG, "est.pdf"), p, width = 9, height = 6)
  cat("  Created est.pdf\n")

  # n1_19.pdf — Small establishments (<20 employees)
  est_did <- feols(n1_19 ~ i(year, treated_clothing, 1977) + i(year, treated, 1977) | unit + year,
                   cluster = ~stcode,
                   data = data[year >= 1974 & year < 1984])
  p <- pretty_3did(est_did, 1977, "Count", "Establishments <20 emp")
  ggsave(file.path(OUTPUT_FIG, "n1_19.pdf"), p, width = 9, height = 6)
  cat("  Created n1_19.pdf\n")

  # frac_small.pdf — Fraction small
  est_did <- feols(frac_small ~ i(year, treated_clothing, 1977) + i(year, treated, 1977) | unit + year,
                   cluster = ~stcode,
                   data = data[year >= 1974 & year < 1984])
  p <- pretty_3did(est_did, 1977, "Fraction", "Fraction Small")
  ggsave(file.path(OUTPUT_FIG, "frac_small.pdf"), p, width = 9, height = 6)
  cat("  Created frac_small.pdf\n")

  # === Simple DiD (no clothing interaction) ===
  est_did <- feols(est ~ i(year, treated, 1977) | unit + year,
                   cluster = ~stcode,
                   data = data[year >= 1974 & year < 1984])
  p <- pretty_did(est_did, 1977, "Count", "Establishments")
  ggsave(file.path(OUTPUT_FIG, "dif_est.pdf"), p, width = 9, height = 6)

  est_did <- feols(n1_19 ~ i(year, treated, 1977) | unit + year,
                   cluster = ~stcode,
                   data = data[year >= 1974 & year < 1984])
  p <- pretty_did(est_did, 1977, "Count", "Establishments <20 emp")
  ggsave(file.path(OUTPUT_FIG, "dif_n1_19.pdf"), p, width = 9, height = 6)

  est_did <- feols(frac_small ~ i(year, treated, 1977) | unit + year,
                   cluster = ~stcode,
                   data = data[year >= 1974 & year < 1984])
  p <- pretty_did(est_did, 1977, "Fraction", "Fraction Small")
  ggsave(file.path(OUTPUT_FIG, "dif_frac_small.pdf"), p, width = 9, height = 6)
  cat("  Created simple DiD figures.\n")

  rm(data); gc()
} else {
  cat("  CBP data not found at:", cbp_file, "\n")
}

cat("BEA/CBP DiD analysis complete.\n")
