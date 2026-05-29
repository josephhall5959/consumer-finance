# 06_clean_bea.R — Clean Bureau of Economic Analysis data

cat("Cleaning BEA data...\n")

# GDP — Real state GDP
gdp_files <- list.files(DATA_RAW, pattern = "SAGDP2S.*\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(gdp_files) == 0) {
  gdp_files <- list.files(file.path(DATA_RAW, "BEA"), pattern = "SAGDP.*\\.csv$", full.names = TRUE)
}

if (length(gdp_files) > 0) {
  data_gdp_list <- lapply(gdp_files, function(x) {
    d <- fread(x, check.names = TRUE)
    # Convert all year columns to character to avoid type mismatch on rbind
    yr_cols <- grep("^X\\d{4}$", names(d), value = TRUE)
    for (col in yr_cols) d[[col]] <- as.character(d[[col]])
    d
  })
  data.in <- rbindlist(data_gdp_list, fill = TRUE)

  # CPI for deflation
  cpi_file <- list.files(file.path(DATA_RAW, "FRED"), pattern = "CPIAUCSL", full.names = TRUE)
  if (length(cpi_file) > 0) {
    cpi <- fread(cpi_file[1]) %>%
      mutate(year = as.integer(substr(DATE, 1, 4)),
             cpi = as.numeric(get(names(.)[2])) / 100) %>%
      dplyr::select(year, cpi) %>%
      group_by(year) %>% summarise(cpi = mean(cpi), .groups = "drop")
  } else {
    cpi <- data.frame(year = 1960:2000, cpi = 1)
  }

  # Find year columns
  year_cols <- grep("^X\\d{4}$", names(data.in), value = TRUE)

  data_gdp <- data.in %>%
    dplyr::select(GeoFIPS, GeoName, Description, all_of(year_cols)) %>%
    pivot_longer(all_of(year_cols), names_to = "Xyear", values_to = "GDP") %>%
    mutate(year = as.numeric(substr(Xyear, 2, 5)),
           fipstate = as.numeric(GeoFIPS) / 1000,
           GDP = suppressWarnings(as.numeric(GDP))) %>%
    filter(!is.na(GDP)) %>%
    dplyr::select(-Xyear) %>%
    left_join(cpi, by = "year") %>%
    mutate(GDP = GDP * cpi)
} else {
  data_gdp <- data.frame()
  warning("No BEA GDP files found")
}

# Employment
emp_files <- list.files(DATA_RAW, pattern = "SAEMP27S.*\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(emp_files) == 0) {
  emp_files <- list.files(file.path(DATA_RAW, "BEA"), pattern = "SAEMP.*\\.csv$", full.names = TRUE)
}

if (length(emp_files) > 0) {
  data.in <- rbindlist(lapply(emp_files, function(x) {
    d <- fread(x, check.names = TRUE)
    yr_cols <- grep("^X\\d{4}$", names(d), value = TRUE)
    for (col in yr_cols) d[[col]] <- as.character(d[[col]])
    d
  }), fill = TRUE)
  year_cols <- grep("^X\\d{4}$", names(data.in), value = TRUE)

  data_emp <- data.in %>%
    dplyr::select(GeoFIPS, GeoName, Description, all_of(year_cols)) %>%
    pivot_longer(all_of(year_cols), names_to = "Xyear", values_to = "Employment") %>%
    mutate(year = as.numeric(substr(Xyear, 2, 5)),
           fipstate = as.numeric(GeoFIPS) / 1000,
           Employment = suppressWarnings(as.numeric(Employment)),
           Description = case_when(
             Description == "Wage and salary employment by place of work (number of jobs)" ~ "All industry total",
             TRUE ~ Description)) %>%
    dplyr::select(-Xyear)
} else {
  data_emp <- data.frame()
}

# Wages
wage_files <- list.files(DATA_RAW, pattern = "SAINC7S.*\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(wage_files) == 0) {
  wage_files <- list.files(file.path(DATA_RAW, "BEA"), pattern = "SAINC.*\\.csv$", full.names = TRUE)
}

if (length(wage_files) > 0) {
  data.in <- rbindlist(lapply(wage_files, function(x) {
    d <- fread(x, check.names = TRUE)
    yr_cols <- grep("^X\\d{4}$", names(d), value = TRUE)
    for (col in yr_cols) d[[col]] <- as.character(d[[col]])
    d
  }), fill = TRUE)
  year_cols <- grep("^X\\d{4}$", names(data.in), value = TRUE)

  data_wage <- data.in %>%
    dplyr::select(GeoFIPS, GeoName, Description, all_of(year_cols)) %>%
    pivot_longer(all_of(year_cols), names_to = "Xyear", values_to = "Wages") %>%
    mutate(year = as.numeric(substr(Xyear, 2, 5)),
           fipstate = as.numeric(GeoFIPS) / 1000,
           Wages = suppressWarnings(as.numeric(Wages)),
           Description = case_when(
             Description == "Wages and salaries by place of work (thousands of dollars)" ~ "All industry total",
             TRUE ~ Description)) %>%
    dplyr::select(-Xyear)
} else {
  data_wage <- data.frame()
}

# State-level prices (Nakamura-Steinsson fiscal stimulus data)
prices_file <- file.path(DATA_RAW, "fiscal_stimulus_coded.dta")
if (file.exists(prices_file)) {
  data_prices <- read_dta(prices_file) %>%
    dplyr::select(year, state, cpi_delneg) %>%
    filter(!is.na(cpi_delneg))
} else {
  data_prices <- data.frame(year = integer(), state = character(), cpi_delneg = numeric())
}

# State IDs
states <- fread(file.path(DATA_RAW, "states.csv"))

# Merge all BEA data
if (nrow(data_gdp) > 0 || nrow(data_emp) > 0 || nrow(data_wage) > 0) {
  data <- data_gdp
  if (nrow(data_emp) > 0) data <- full_join(data, data_emp, by = c("GeoFIPS", "GeoName", "Description", "year", "fipstate"))
  if (nrow(data_wage) > 0) data <- full_join(data, data_wage, by = c("GeoFIPS", "GeoName", "Description", "year", "fipstate"))

  # Ensure columns exist even if some data sources are empty
  if (!"GDP" %in% names(data)) data$GDP <- NA_real_
  if (!"Employment" %in% names(data)) data$Employment <- NA_real_
  if (!"Wages" %in% names(data)) data$Wages <- NA_real_

  data <- data %>%
    mutate(across(c(any_of(c("Wages", "Employment", "GDP"))), as.numeric)) %>%
    distinct() %>%
    left_join(states %>% dplyr::select(fipstate = st, state = Code), by = "fipstate") %>%
    left_join(data_prices, by = c("state", "year")) %>%
    group_by(fipstate, Description) %>%
    arrange(fipstate, Description, year) %>%
    mutate(
      GDP = ifelse(!is.na(cpi_delneg) & cpi_delneg != 0, GDP / cpi_delneg, GDP),
      Wages = ifelse(!is.na(cpi_delneg) & cpi_delneg != 0, Wages / cpi_delneg, Wages),
      wage_rate = Wages / Employment,
      gdp_growth = GDP / lag(GDP) - 1,
      emp_growth = Employment / lag(Employment) - 1,
      wage_growth = wage_rate / lag(wage_rate) - 1,
      price_growth = cpi_delneg / lag(cpi_delneg) - 1
    ) %>%
    filter(year > 1969 & year < 1991 & fipstate > 0 & fipstate < 90) %>%
    ungroup() %>%
    as.data.table()

  # Industry categories
  industry_cat_file <- file.path(DATA_RAW, "BEA", "industry_cat.xlsx")
  if (file.exists(industry_cat_file)) {
    industry_cat <- read_excel(industry_cat_file)
    data[, Description := trimws(Description)]
    data <- merge(data, industry_cat, by = "Description", all.x = TRUE)
  }

  # Add usury treatment
  usury <- fread(file.path(DATA_GEN, "usury.csv")) %>%
    dplyr::select(fipstate = st, r_77) %>%
    distinct() %>%
    mutate(treated = ifelse(r_77 < 18, 1, 0))

  data <- left_join(data, usury, by = "fipstate")

  fwrite(data, file.path(DATA_GEN, "bea.csv"))
  cat("BEA cleaning complete.\n")
} else {
  warning("No BEA data found to process")
  fwrite(data.table(), file.path(DATA_GEN, "bea.csv"))
}
