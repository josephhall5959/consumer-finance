# 05_clean_psid.R — Clean Panel Study of Income Dynamics data

cat("Cleaning PSID data...\n")

psid_excel <- file.path(DATA_RAW, "PSID", "J308439.xlsx")
psid_available <- file.exists(psid_excel) && file.size(psid_excel) > 0

if (!psid_available) {
  # Try alternative filenames (non-zero only)
  psid_files <- list.files(file.path(DATA_RAW, "PSID"), pattern = "\\.xlsx$", full.names = TRUE)
  psid_files <- psid_files[file.size(psid_files) > 0]
  if (length(psid_files) > 0) {
    psid_excel <- psid_files[1]
    psid_available <- TRUE
  }
}

if (!psid_available) {
  stop("No PSID data found. Place J308439.xlsx in ", file.path(DATA_RAW, "PSID"))
} else {
  data.in <- read_excel(psid_excel)

  codebook_file <- file.path(DATA_RAW, "PSID", "psid.xlsx")
  varlabels_file <- file.path(DATA_RAW, "PSID", "varlabels_in.xlsx")

  if (!file.exists(varlabels_file) || file.size(varlabels_file) == 0) {
    warning("PSID varlabels_in.xlsx not found. Skipping PSID cleaning.")
    fwrite(data.table(), file.path(DATA_GEN, "psid_final.csv"))
    cat("PSID cleaning skipped (no varlabels).\n")
  } else {
    codebook <- read_excel(codebook_file)
    varlabels <- read_excel(varlabels_file)

    index <- codebook %>%
      mutate(series_number = 1:n()) %>%
      dplyr::select(series_number, starts_with("Y")) %>%
      pivot_longer(!series_number, names_to = "year", names_prefix = "Y")

    data.long <- data.in %>%
      mutate(id = 1:n()) %>%
      pivot_longer(-id, names_to = "name") %>%
      left_join(varlabels, by = "name")

    data.wide <- data.long %>%
      filter(varname != "family_number") %>%
      dplyr::select(-name, -series_number) %>%
      pivot_wider(names_from = varname, values_from = value) %>%
      filter(!is.na(state_psid))

    # Merge state IDs
    states <- fread(file.path(DATA_RAW, "states.csv")) %>%
      dplyr::select(state_psid, State = Code)

    data <- data.wide %>% left_join(states, by = "state_psid")

    # Attach usury data
    usury_rates <- fread(file.path(DATA_RAW, "usury_rates.csv")) %>%
      dplyr::select(State, y1977)

    data <- data %>%
      left_join(usury_rates, by = "State") %>%
      as.data.table()

    # Construct key variables
    data[, food_insufficient := ifelse(food_expense < food_req, 1, 0)]
    data[, black := ifelse(race == 2, 1, 0)]
    data[, hispanic := ifelse(race == 3, 1, 0)]
    data[, nonwhite := ifelse(race != 1, 1, 0)]
    data[, owns := case_when(homeowner == 1 ~ 1, TRUE ~ 0)]
    data[, rents := case_when(homeowner == 5 ~ 1, TRUE ~ 0)]
    data[, move := ifelse(lead(State) != State, 1, 0), by = id]
    data[, change_hh_size := hh_size - lag(hh_size), by = id]
    data[, ncars := ifelse(ncars > 8, NA_real_, ncars)]
    data[, retirement := ifelse(retired == 1 & lag(retired) == 0, 1, 0), by = id]
    data[, treated := case_when(y1977 %in% c("12", "15") ~ 1, TRUE ~ 0)]

    # Clean bankruptcy
    bk1 <- data %>%
      group_by(id, bk1_year) %>% summarise(n = n(), .groups = "drop") %>%
      filter(bk1_year > 0 & bk1_year < 2000) %>%
      dplyr::select(id, year = bk1_year) %>%
      mutate(bk1 = 1)

    bk2 <- data %>%
      group_by(id, bk2_year) %>% summarise(n = n(), .groups = "drop") %>%
      filter(bk2_year > 0 & bk2_year < 2000) %>%
      dplyr::select(id, year = bk2_year) %>%
      mutate(bk2 = 1)

    data <- left_join(data, bk1, by = c("id", "year")) %>%
      left_join(bk2, by = c("id", "year")) %>%
      mutate(bankrupt = case_when(bk1 == 1 ~ 1, bk2 == 1 ~ 1, TRUE ~ 0))

    # Clean outliers
    data[, income := ifelse(income < 9999999, income, NA_real_)]
    data[, childcare_expense := ifelse(childcare_expense < 99999, childcare_expense, NA_real_)]

    fwrite(data, file.path(DATA_GEN, "psid_final.csv"))
    cat("PSID cleaning complete.\n")
  }
}
