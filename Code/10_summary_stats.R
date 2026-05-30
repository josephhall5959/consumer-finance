# 10_summary_stats.R — Summary statistics and time-series figures

cat("Generating summary statistics...\n")

# Load cleaned data
scf_file <- file.path(DATA_GEN, "clean_ts.csv")
if (!file.exists(scf_file) || file.info(scf_file)$size < 100) {
  warning("clean_ts.csv not found. Some summary stats will be skipped.")
  scf <- data.table()
} else {
  scf <- fread(scf_file)
}

usury_file <- file.path(DATA_GEN, "usury.csv")
if (file.exists(usury_file)) {
  usury <- fread(usury_file)
} else {
  usury <- fread(file.path(DATA_RAW, "usury_rates.csv"))
}

# --- SCF Summary Statistics ---
if (nrow(scf) > 0) {

  # Card ownership over time
  card_ts <- scf %>%
    group_by(year) %>%
    summarise(
      has_card = weighted.mean(has_card, weight, na.rm = TRUE),
      has_bank_card = weighted.mean(has_bank_card, weight, na.rm = TRUE),
      has_store_card = weighted.mean(has_store_card, weight, na.rm = TRUE),
      has_gas_card = weighted.mean(ifelse(is.na(has_gas_card), 0, has_gas_card), weight, na.rm = TRUE),
      card_spend_mean = weighted.mean(card_spend, weight, na.rm = TRUE),
      card_bal_mean = weighted.mean(card_bal, weight, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # Figure: Spending levels (fraction of income)
  if ("card_spend" %in% names(scf)) {
    yearly_summary <- scf[, .(
      store_card_spend = weighted.mean(store_card_spend, w = weight, na.rm = TRUE),
      bank_card_spend = weighted.mean(bank_card_spend, w = weight, na.rm = TRUE),
      income = weighted.mean(income, w = weight, na.rm = TRUE)
    ), by = "year"]
    yearly_summary[, store_to_income := 12 * store_card_spend / income]
    yearly_summary[, bank_to_income := 12 * bank_card_spend / income]

    toplot <- yearly_summary[!year %in% c(1977, 1983, 1986)] %>%
      dplyr::select(year, `Store Card Spending` = store_to_income,
                     `Bank Card Spending` = bank_to_income) %>%
      pivot_longer(-year, names_to = "Series", values_to = "value")

    p <- ggplot(toplot, aes(x = year, y = value, color = Series, linetype = Series)) +
      geom_line() +
      labs(x = "Year", y = "Fraction of Income") +
      scale_y_continuous(labels = scales::percent, limits = c(0, .2)) +
      theme_classic() +
      theme(legend.position = c(0.2, 0.8))
    ggsave(file.path(OUTPUT_FIG, "spending_levels.pdf"), p, width = 9, height = 6)
  }

  # Figure: Store card fraction
  if ("has_store_card" %in% names(scf) && "has_bank_card" %in% names(scf)) {
    frac_ts <- card_ts %>%
      mutate(store_frac = has_store_card / (has_store_card + has_bank_card))

    p <- ggplot(frac_ts, aes(x = year, y = store_frac)) +
      geom_line() + geom_point() +
      theme_classic() +
      labs(x = "Year", y = "Store card share of cardholders") +
      scale_y_continuous(limits = c(0, 1))
    ggsave(file.path(OUTPUT_FIG, "store_fraction.pdf"), p, width = 9, height = 6)
  }

  # Figure: Cards per household (bank vs store), SCF
  # Restrict to 1977 onward (the Marquette era) and drop the 1983 survey, whose
  # pre-modern-SCF sampling produces an erratic reading inconsistent with the
  # methodologically-comparable triennial SCF that begins in 1989.
  if (all(c("bank_cards", "store_cards") %in% names(scf))) {
    cards_ts <- scf[, .(
      `Bank Cards`  = weighted.mean(bank_cards,  w = weight, na.rm = TRUE),
      `Store Cards` = weighted.mean(store_cards, w = weight, na.rm = TRUE)
    ), by = "year"]
    cards_ts <- cards_ts[!is.na(year) & year >= 1977 & year != 1983 &
                           !is.nan(`Bank Cards`) & !is.nan(`Store Cards`)]
    cards_long <- cards_ts %>%
      pivot_longer(-year, names_to = "Series", values_to = "value")
    p <- ggplot(cards_long, aes(x = year, y = value, color = Series, linetype = Series)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.3) +
      scale_color_manual(values = c("Bank Cards" = "firebrick", "Store Cards" = "steelblue")) +
      labs(x = "Year", y = "Cards per household") +
      theme_classic() +
      theme(legend.position = c(0.25, 0.85), legend.title = element_blank())
    ggsave(file.path(OUTPUT_FIG, "cards_per_household.pdf"), p, width = 9, height = 6)
  }

  # Figure: Bank card ownership
  bank_card_ts <- scf[, .(has_bank_card = mean(has_bank_card)), by = year] %>%
    rbind(data.frame(year = c(1968, 2020), has_bank_card = c(0.07, 0.8)))
  p <- ggplot(bank_card_ts[!is.na(has_bank_card)], aes(x = year, y = has_bank_card)) +
    geom_line() +
    labs(x = "Year") +
    scale_y_continuous(name = "Households with bank credit card", labels = scales::percent) +
    theme_classic()
  ggsave(file.path(OUTPUT_FIG, "has_bank_card.pdf"), p, width = 9, height = 6)

  # Figure: Card growth stacked (call report data by bank HQ state)
  call_report_file <- file.path(DATA_RAW, "Call_Reports", "call_reports.csv")
  if (!file.exists(call_report_file)) {
    call_report_file <- file.path(DATA_GEN, "call_reports.csv")
  }
  if (file.exists(call_report_file) && file.info(call_report_file)$size > 100) {
    cr <- fread(call_report_file) %>% as_tibble()
    all_cards <- cr %>% group_by(year) %>%
      summarise(cards = sum(cards, na.rm = TRUE) / 1e6, .groups = "drop")
    dereg_cards <- cr %>%
      filter(state %in% c("DE", "SD", "AZ", "IL", "ID", "MT", "NV",
                           "NH", "NJ", "NM", "OR", "UT", "WI")) %>%
      group_by(year) %>%
      summarise(cards_dereg = sum(cards, na.rm = TRUE) / 1e6, .groups = "drop")
    toplot_cr <- merge(all_cards, dereg_cards) %>%
      mutate(cards_reg = cards - cards_dereg)

    gdp_file <- file.path(DATA_RAW, "GDP.csv")
    if (file.exists(gdp_file)) {
      gdp <- fread(gdp_file) %>%
        mutate(year = as.integer(substr(DATE, 1, 4)), gdp = as.numeric(GDP))
      toplot_pct <- toplot_cr %>% merge(gdp) %>%
        mutate(Regulated = cards_reg / gdp, Deregulated = cards_dereg / gdp) %>%
        dplyr::select(year, Regulated, Deregulated) %>%
        pivot_longer(-year, names_to = "Bank HQ State", values_to = "Card Lending/GDP")

      p <- ggplot(toplot_pct, aes(x = year, y = `Card Lending/GDP`, fill = `Bank HQ State`)) +
        geom_bar(stat = "identity") +
        theme_classic()
      ggsave(file.path(OUTPUT_FIG, "card_growth_stacked.pdf"), p, width = 6, height = 4)
    }
    rm(cr, all_cards, dereg_cards, toplot_cr); gc()
  }

  # Figure: Interest rate distribution (rates_all.png)
  if ("card_rate" %in% names(scf)) {
    rate_data <- scf %>% filter(!is.na(card_rate) & card_rate > 0 & card_rate < 50)
    if (nrow(rate_data) > 0) {
      p <- ggplot(rate_data, aes(x = card_rate)) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
        theme_classic() +
        labs(x = "Interest rate", y = "Count") +
        facet_wrap(~year, scales = "free_y")
      ggsave(file.path(OUTPUT_FIG, "rates_all.png"), p, width = 12, height = 8, dpi = 300)
    }
  }

  # Figure: Bunching at usury limits (raw card rate histograms)
  if ("card_rate" %in% names(scf)) {
    bunching_77 <- scf %>% filter(year == 1977 & card_rate < 40)
    bunching_83 <- scf %>% filter(year == 1983 & card_rate < 40)

    if (nrow(bunching_77) > 0) {
      p <- ggplot(bunching_77, aes(x = card_rate)) + geom_histogram() +
        labs(x = "1977 Card Rate") + theme_classic()
      ggsave(file.path(OUTPUT_FIG, "bunching.pdf"), p, width = 6, height = 6)
    }
    if (nrow(bunching_83) > 0) {
      p <- ggplot(bunching_83, aes(x = card_rate)) + geom_histogram() +
        labs(x = "1983 Card Rate") + theme_classic()
      ggsave(file.path(OUTPUT_FIG, "bunching_83.pdf"), p, width = 6, height = 6)
    }
  }

  # Table: SCF summary statistics (summary_scf.tex)
  summary_vars <- c("income", "has_card", "has_bank_card", "has_store_card",
                     "card_spend", "card_bal")
  avail_vars <- intersect(summary_vars, names(scf))

  if (length(avail_vars) > 0) {
    sum_tab <- scf %>%
      dplyr::select(all_of(avail_vars)) %>%
      summarise(across(everything(), list(
        Mean = ~mean(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE),
        N = ~sum(!is.na(.x))
      ))) %>%
      pivot_longer(everything(),
                   names_to = c("Variable", ".value"),
                   names_pattern = "(.+)_(.+)")

    sum_tab %>%
      kbl(format = "latex", booktabs = TRUE, digits = 2,
          caption = "Summary Statistics: Survey of Consumer Finances") %>%
      kable_classic() %>%
      save_kable(file = file.path(OUTPUT_TAB, "summary_scf.tex"),
                 self_contained = FALSE)
  }

  # Table: Card uses (card_uses.tex) - from comparison 1970.do
  if (all(c("has_bank_card", "has_store_card", "year") %in% names(scf))) {
    card_uses <- scf %>%
      filter(year %in% c(1970, 1977, 1983, 2019)) %>%
      group_by(year) %>%
      summarise(
        `Bank card` = weighted.mean(has_bank_card, weight, na.rm = TRUE),
        `Store card` = weighted.mean(has_store_card, weight, na.rm = TRUE),
        `Any card` = weighted.mean(has_card, weight, na.rm = TRUE),
        N = n(),
        .groups = "drop"
      )

    # Bare tabular only: main.tex inputs this inside a \resizebox under a
    # shared caption, so an inner table environment would break compilation.
    cu_tex <- card_uses %>%
      kbl(format = "latex", booktabs = TRUE, digits = 2) %>%
      kable_classic() %>%
      as.character()
    cu_lines <- strsplit(cu_tex, "\n")[[1]]
    cu_lines <- cu_lines[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption", cu_lines)]
    writeLines(cu_lines, file.path(OUTPUT_TAB, "card_uses.tex"))
  }
}

# --- Inflation time series ---
fedfunds_file <- file.path(DATA_RAW, "FEDFUNDS.csv")
inflation_file <- file.path(DATA_RAW, "FPCPITOTLZGUSA.csv")

if (file.exists(fedfunds_file) && file.exists(inflation_file)) {
  ff <- fread(fedfunds_file) %>%
    mutate(year = as.integer(substr(DATE, 1, 4)),
           FEDFUNDS = as.numeric(FEDFUNDS)) %>%
    group_by(year) %>%
    summarise(fedfunds = mean(FEDFUNDS, na.rm = TRUE), .groups = "drop")

  inf <- fread(inflation_file)
  inf_val_col <- setdiff(names(inf), "DATE")[1]
  inf <- inf %>%
    mutate(year = as.integer(substr(DATE, 1, 4)),
           inflation = as.numeric(get(inf_val_col))) %>%
    dplyr::select(year, inflation)

  macro_ts <- full_join(ff, inf, by = "year") %>% filter(year >= 1960 & year <= 2000)

  macro_long <- macro_ts %>%
    pivot_longer(-year, names_to = "series", values_to = "value")

  p <- ggplot(macro_long, aes(x = year, y = value, color = series)) +
    geom_line() +
    theme_classic() +
    labs(x = "Year", y = "Percent", color = "") +
    theme(legend.position = "bottom")
  ggsave(file.path(OUTPUT_FIG, "inflation_ts.pdf"), p, width = 9, height = 6)
}

# --- International comparison (scatter plot) ---
intl_file <- file.path(DATA_RAW, "cross_country.xlsx")
if (file.exists(intl_file)) {
  intl <- read_excel(intl_file, sheet = "Sheet1")

  # Find the label and key columns
  label_col <- grep("Label", names(intl), value = TRUE)[1]
  cc_col <- grep("Credit.Card.Ownership|Credit Card Ownership", names(intl), value = TRUE)[1]
  sc_col <- grep("Store.Credit.Borrowing|Store Credit Borrowing", names(intl), value = TRUE)[1]

  if (!is.na(cc_col) && !is.na(sc_col)) {
    # Scatter plot: use only obs with both variables
    toplot_intl <- intl %>%
      filter(!is.na(.data[[cc_col]]) & !is.na(.data[[sc_col]]))

    lbl <- if (!is.na(label_col)) toplot_intl[[label_col]] else ""

    p <- ggplot(toplot_intl, aes(x = .data[[cc_col]], y = .data[[sc_col]])) +
      geom_point() +
      geom_text(label = lbl, nudge_y = -0.013, nudge_x = 0.04) +
      theme_classic()
    ggsave(file.path(OUTPUT_FIG, "international.pdf"), p, width = 9, height = 6)
    cat("  Created international.pdf\n")

    # Histograms (for slides): keep latest year per country, require both vars
    country_col <- grep("Country", names(intl), value = TRUE)[1]
    year_col <- grep("Year", names(intl), value = TRUE)[1]
    hist_data <- intl %>%
      filter(!is.na(.data[[cc_col]]) & !is.na(.data[[sc_col]])) %>%
      arrange(.data[[country_col]], desc(.data[[year_col]])) %>%
      group_by(.data[[country_col]]) %>%
      slice(1) %>%
      ungroup()

    # US and Canada values for vertical lines
    us_cc <- hist_data %>% filter(.data[[country_col]] == "United States") %>% pull(!!cc_col)
    us_sc <- hist_data %>% filter(.data[[country_col]] == "United States") %>% pull(!!sc_col)
    ca_cc <- hist_data %>% filter(.data[[country_col]] == "Canada") %>% pull(!!cc_col)
    ca_sc <- hist_data %>% filter(.data[[country_col]] == "Canada") %>% pull(!!sc_col)

    p1 <- ggplot(hist_data, aes(x = .data[[cc_col]])) +
      geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = us_cc, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = us_cc, y = Inf, label = "US", vjust = 2, hjust = -0.3, color = "red") +
      theme_classic() +
      labs(x = "Credit Card Ownership", y = "Number of Countries", title = "(a) Credit Card Ownership")
    if (length(ca_cc) > 0 && !is.na(ca_cc)) {
      p1 <- p1 +
        geom_vline(xintercept = ca_cc, linetype = "dashed", color = "darkgreen", linewidth = 1) +
        annotate("text", x = ca_cc, y = Inf, label = "Canada", vjust = 3.5, hjust = 1.15, color = "darkgreen")
    }

    p2 <- ggplot(hist_data, aes(x = .data[[sc_col]])) +
      geom_histogram(bins = 20, fill = "coral", alpha = 0.7) +
      geom_vline(xintercept = us_sc, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = us_sc, y = Inf, label = "US", vjust = 2, hjust = -0.3, color = "red") +
      theme_classic() +
      labs(x = "Store Credit Borrowing", y = "Number of Countries", title = "(b) Store Credit Borrowing")
    if (length(ca_sc) > 0 && !is.na(ca_sc)) {
      p2 <- p2 +
        geom_vline(xintercept = ca_sc, linetype = "dashed", color = "darkgreen", linewidth = 1) +
        annotate("text", x = ca_sc, y = Inf, label = "Canada", vjust = 2, hjust = 1.15, color = "darkgreen")
    }

    p_combined <- grid.arrange(p1, p2, ncol = 2)
    ggsave(file.path(OUTPUT_FIG, "international_histograms.pdf"), p_combined, width = 12, height = 5)
    cat(sprintf("  Created international_histograms.pdf (%d countries)\n", nrow(hist_data)))
  }
}

# Clean up SCF before loading CBP
rm(scf); gc()

# --- CBP summary statistics ---
cbp_file <- file.path(DATA_GEN, "cbp_sic2.csv")
if (file.exists(cbp_file) && file.info(cbp_file)$size > 100) {
  cbp <- fread(cbp_file, select = c("year", "sic2", "emp", "est", "n1_4", "n1_19", "state"))

  summary_tab <- cbp %>%
    summarise(across(where(is.numeric), list(
      Mean = ~mean(.x, na.rm = TRUE),
      SD = ~sd(.x, na.rm = TRUE),
      N = ~sum(!is.na(.x))
    ))) %>%
    pivot_longer(everything(),
                 names_to = c("Variable", ".value"),
                 names_pattern = "(.+)_(.+)")

  summary_tab %>%
    kbl(format = "latex", booktabs = TRUE, digits = 2,
        caption = "Summary Statistics: County Business Patterns") %>%
    kable_classic() %>%
    save_kable(file = file.path(OUTPUT_TAB, "summary_cbp.tex"),
               self_contained = FALSE)
  rm(cbp); gc()
}

# --- A5: Treatment balance table ---
# Compare pre-treatment (1977) state characteristics between treated and control states
usury_file2 <- file.path(DATA_GEN, "usury.csv")
bea_file2 <- file.path(DATA_GEN, "bea.csv")
cbp_file2 <- file.path(DATA_GEN, "cbp_sic2.csv")

if (file.exists(usury_file2) && file.exists(bea_file2)) {
  usury_bal <- fread(usury_file2) %>%
    dplyr::select(state, r_77) %>% distinct() %>%
    mutate(treated = ifelse(r_77 < 18, 1, 0))

  bea_bal <- fread(bea_file2)
  bea_77 <- bea_bal %>%
    filter(year == 1977) %>%
    group_by(state) %>%
    summarise(
      gdp_total = sum(GDP, na.rm = TRUE),
      emp_total = sum(Employment, na.rm = TRUE),
      wages_total = sum(Wages, na.rm = TRUE),
      wage_rate = mean(wage_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(usury_bal, by = "state") %>%
    filter(!is.na(treated))

  # If CBP available, add retail employment share
  if (file.exists(cbp_file2) && file.info(cbp_file2)$size > 100) {
    cbp_bal <- fread(cbp_file2, select = c("year", "state", "sic2", "emp"))
    cbp_77 <- cbp_bal %>%
      filter(year == 1977) %>%
      group_by(state) %>%
      summarise(total_emp = sum(emp, na.rm = TRUE), .groups = "drop")
    cbp_retail_77 <- cbp_bal %>%
      filter(year == 1977 & sic2 >= 52 & sic2 <= 59) %>%
      group_by(state) %>%
      summarise(retail_emp = sum(emp, na.rm = TRUE), .groups = "drop")
    cbp_share <- cbp_77 %>%
      left_join(cbp_retail_77, by = "state") %>%
      mutate(retail_share = retail_emp / total_emp)
    bea_77 <- bea_77 %>% left_join(cbp_share %>% dplyr::select(state, retail_share), by = "state")
    rm(cbp_bal, cbp_77, cbp_retail_77, cbp_share); gc()
  }

  # Build balance table
  balance_vars <- intersect(c("wage_rate", "emp_total", "gdp_total", "retail_share"), names(bea_77))
  balance_labels <- c(wage_rate = "Wage rate", emp_total = "Employment",
                      gdp_total = "GDP", retail_share = "Retail emp. share")

  balance_rows <- lapply(balance_vars, function(v) {
    t_vals <- bea_77 %>% filter(treated == 1) %>% pull(!!sym(v))
    c_vals <- bea_77 %>% filter(treated == 0) %>% pull(!!sym(v))
    tt <- tryCatch(t.test(t_vals, c_vals), error = function(e) NULL)
    data.frame(
      Variable = balance_labels[v],
      `Treated Mean` = round(mean(t_vals, na.rm = TRUE), 2),
      `Control Mean` = round(mean(c_vals, na.rm = TRUE), 2),
      Difference = round(mean(t_vals, na.rm = TRUE) - mean(c_vals, na.rm = TRUE), 2),
      `p-value` = if (!is.null(tt)) round(tt$p.value, 3) else NA,
      check.names = FALSE
    )
  })
  balance_tab <- do.call(rbind, balance_rows)

  tex_out <- balance_tab %>%
    kbl(format = "latex", booktabs = TRUE) %>%
    kable_classic() %>%
    as.character()
  tex_lines <- strsplit(tex_out, "\n")[[1]]
  tex_lines <- tex_lines[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption", tex_lines)]
  writeLines(tex_lines, file.path(OUTPUT_TAB, "treatment_balance.tex"))
  cat("  Created treatment_balance.tex\n")
  rm(bea_bal, bea_77, usury_bal); gc()
}

# --- A12: Card rates raw series ---
# Show raw bank vs store card interest rate series by year
scf_file_a12 <- file.path(DATA_GEN, "clean_ts.csv")
if (file.exists(scf_file_a12) && file.info(scf_file_a12)$size > 100) {
  scf_a12 <- fread(scf_file_a12)
} else {
  scf_a12 <- data.table()
}

if (nrow(scf_a12) > 0 && "card_rate" %in% names(scf_a12)) {
  # Attempt to distinguish bank vs store card rates
  rate_cols <- grep("rate", names(scf_a12), value = TRUE, ignore.case = TRUE)
  bank_rate_col <- grep("bank.*rate|rate.*bank", names(scf_a12), value = TRUE, ignore.case = TRUE)[1]
  store_rate_col <- grep("store.*rate|rate.*store", names(scf_a12), value = TRUE, ignore.case = TRUE)[1]

  if (!is.na(bank_rate_col) && !is.na(store_rate_col)) {
    rate_ts <- scf_a12 %>%
      group_by(year) %>%
      summarise(
        bank_rate = weighted.mean(get(bank_rate_col), weight, na.rm = TRUE),
        store_rate = weighted.mean(get(store_rate_col), weight, na.rm = TRUE),
        bank_zero_frac = weighted.mean(get(bank_rate_col) == 0, weight, na.rm = TRUE),
        store_zero_frac = weighted.mean(get(store_rate_col) == 0, weight, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    # Fallback: use card_rate overall and zero-rate fraction
    rate_ts <- scf_a12 %>%
      filter(!is.na(card_rate)) %>%
      group_by(year) %>%
      summarise(
        mean_rate = weighted.mean(card_rate, weight, na.rm = TRUE),
        zero_frac = weighted.mean(card_rate == 0, weight, na.rm = TRUE),
        .groups = "drop"
      )
  }

  if ("bank_rate" %in% names(rate_ts)) {
    rate_long <- rate_ts %>%
      dplyr::select(year, `Bank card` = bank_rate, `Store card` = store_rate) %>%
      pivot_longer(-year, names_to = "type", values_to = "rate")
    zero_long <- rate_ts %>%
      dplyr::select(year, `Bank card` = bank_zero_frac, `Store card` = store_zero_frac) %>%
      pivot_longer(-year, names_to = "type", values_to = "zero_frac")

    p1 <- ggplot(rate_long, aes(x = year, y = rate, color = type)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      theme_classic(base_size = 14) +
      labs(x = "Year", y = "Mean interest rate (%)", color = "Card type",
           subtitle = "(a) Mean interest rate") +
      theme(legend.position = "bottom", axis.title = element_text(size = 13))

    p2 <- ggplot(zero_long, aes(x = year, y = zero_frac, color = type)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      theme_classic(base_size = 14) +
      labs(x = "Year", y = "Fraction with zero rate", color = "Card type",
           subtitle = "(b) Zero-rate fraction") +
      theme(legend.position = "bottom", axis.title = element_text(size = 13))

    p_rates <- gridExtra::grid.arrange(p1, p2, ncol = 2)
    ggsave(file.path(OUTPUT_FIG, "card_rates_raw.pdf"), p_rates, width = 12, height = 6)
    cat("  Created card_rates_raw.pdf\n")
  } else if ("mean_rate" %in% names(rate_ts)) {
    p <- ggplot(rate_ts, aes(x = year)) +
      geom_line(aes(y = mean_rate), linewidth = 1, color = "steelblue") +
      geom_point(aes(y = mean_rate), size = 2, color = "steelblue") +
      geom_line(aes(y = zero_frac * max(rate_ts$mean_rate, na.rm = TRUE)),
                linewidth = 1, color = "firebrick", linetype = "dashed") +
      scale_y_continuous(
        name = "Mean interest rate (%)",
        sec.axis = sec_axis(~ . / max(rate_ts$mean_rate, na.rm = TRUE),
                            name = "Fraction with zero rate")
      ) +
      theme_classic(base_size = 14) +
      labs(x = "Year") +
      theme(axis.title = element_text(size = 13))
    ggsave(file.path(OUTPUT_FIG, "card_rates_raw.pdf"), p, width = 9, height = 6)
    cat("  Created card_rates_raw.pdf\n")
  }
  if (exists("scf_a12")) { rm(scf_a12); gc() }
}

# --- A13: Mian-Sufi-Werner deregulation correlation ---
ks_file <- file.path(DATA_RAW, "kroszner_strahan.xlsx")
if (file.exists(ks_file) && file.exists(usury_file)) {
  ks <- read_excel(ks_file)
  names(ks) <- make.names(names(ks))

  usury_msw <- fread(usury_file) %>%
    dplyr::select(state, r_77) %>% distinct() %>%
    mutate(treated = ifelse(r_77 < 18, 1, 0))

  # Compute MSW-style deregulation index: average year of deregulation
  if (all(c("ma_branch", "full_branch", "interstate") %in% names(ks))) {
    ks <- ks %>%
      mutate(msw_index = (ma_branch + full_branch + interstate) / 3)
  } else {
    # Use whatever columns are available
    dereg_cols <- setdiff(names(ks), "state")
    ks <- ks %>%
      mutate(msw_index = rowMeans(across(all_of(dereg_cols)), na.rm = TRUE))
  }

  msw_merged <- ks %>%
    left_join(usury_msw, by = "state") %>%
    filter(!is.na(r_77) & !is.na(msw_index))

  corr_val <- cor(msw_merged$r_77, msw_merged$msw_index, use = "complete.obs")

  p <- ggplot(msw_merged, aes(x = r_77, y = msw_index)) +
    geom_point(aes(color = factor(treated)), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = "dashed") +
    theme_classic(base_size = 14) +
    labs(x = "1977 usury rate (%)",
         y = "Avg. deregulation year (MSW index)",
         color = "Treated") +
    scale_color_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                       labels = c("Control", "Treated")) +
    annotate("text", x = min(msw_merged$r_77, na.rm = TRUE),
             y = max(msw_merged$msw_index, na.rm = TRUE),
             label = paste0("Corr = ", round(corr_val, 2)),
             hjust = 0, vjust = 1, size = 4.5) +
    theme(legend.position = "bottom", axis.title = element_text(size = 13))
  ggsave(file.path(OUTPUT_FIG, "msw_correlation.pdf"), p, width = 8, height = 6)
  cat("  Created msw_correlation.pdf (correlation =", round(corr_val, 2), ")\n")
  rm(ks, usury_msw, msw_merged); gc()
}

cat("Summary statistics complete.\n")
