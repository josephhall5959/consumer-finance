# 04_clean_cps.R — Clean Current Population Survey data
# Source: IPUMS CPS extract

cat("Cleaning CPS data...\n")

# Check for ipumsr
if (!requireNamespace("ipumsr", quietly = TRUE)) {
  install.packages("ipumsr")
}
library(ipumsr)

states <- fread(file.path(DATA_RAW, "states.csv")) %>% as_tibble()
usury <- fread(file.path(DATA_GEN, "usury.csv")) %>% as_tibble()

# Read IPUMS CPS extract — prefer cps_00010 (has OWNERSHP, PUBHOUS, etc.)
ddi_file <- list.files(file.path(DATA_RAW, "CPS"), pattern = "\\.xml$", full.names = TRUE)
# Filter out 0-byte stub files
ddi_file <- ddi_file[file.size(ddi_file) > 0]
# Sort descending so newer/more complete extracts come first
ddi_file <- sort(ddi_file, decreasing = TRUE)
# But prefer cps_00010 specifically if it exists (has all required variables)
preferred <- grep("cps_00010", ddi_file, value = TRUE)
if (length(preferred) > 0) ddi_file <- c(preferred, setdiff(ddi_file, preferred))

cps_done <- FALSE

if (length(ddi_file) > 0) {
  ddi <- read_ipums_ddi(ddi_file[1])
  data <- read_ipums_micro(ddi)

  data.cps1 <- data %>%
    filter(AGE <= 65 & AGE >= 18 &
             CLASSWKR %!in% c(0, 26, 99) &
             RACE %in% c(100, 200) &
             IND %!in% c(17, 18, 10, 11, 12))

  data.cps2 <- data.cps1 %>%
    mutate(
      black = ifelse(RACE == 200, 1, 0),
      self.employed = ifelse(CLASSWKR %in% c(10, 13, 14), 1, 0),
      ue = ifelse(EMPSTAT != 10, 1, 0),
      emp = ifelse(self.employed == 0 & ue == 0, 1, 0),
      health.shock = ifelse(WHYABSNT %in% c(6, 11) | WHYPTLWK == 100, 1, 0),
      fulltime = ifelse(AHRSWORKT >= 35 & AHRSWORKT != 999, 1, 0),
      rents = ifelse(OWNERSHP == 10, 0, 1),
      poverty = ifelse(OFFPOV == 1, 1, 0),
      married = ifelse(MARST == 1, 1, 0),
      divorced_sep = ifelse(MARST %in% c(3, 4), 1, 0),
      moved = case_when(MIGRATE1 > 1 ~ 1, TRUE ~ 0),
      pub_sub = ifelse(PUBHOUS == 2 | RENTSUB == 2, 1, 0),
      moved.same = case_when(MIGRATE1 == 3 ~ 1, TRUE ~ 0),
      moved.new = case_when(MIGRATE1 %in% c(4, 5) ~ 1, TRUE ~ 0),
      educ.yrs = case_when(
        EDUC %in% c(1, 2) ~ 0, EDUC == 11 ~ 1, EDUC == 12 ~ 2,
        EDUC == 13 ~ 3, EDUC == 14 ~ 4, EDUC == 21 ~ 5, EDUC == 22 ~ 6,
        EDUC == 31 ~ 7, EDUC == 32 ~ 8, EDUC == 40 ~ 9, EDUC == 50 ~ 10,
        EDUC == 60 ~ 11, EDUC %in% c(70, 71, 72, 73) ~ 12, EDUC == 80 ~ 13,
        EDUC == 90 ~ 14, EDUC == 100 ~ 15, EDUC == 110 ~ 16,
        EDUC == 121 ~ 17, EDUC == 122 ~ 18
      )
    ) %>%
    dplyr::select(year = YEAR, id = CPSIDP, fipstate = STATEFIP, fipmsa = METFIPS,
           income = HHINCOME, wt = ASECWT, age = AGE, sex = SEX, everything()) %>%
    left_join(states %>% dplyr::select(state_name = State, state = Code, fipstate = st),
              by = "fipstate") %>%
    left_join(usury, by = c("state", "year")) %>%
    mutate(id = ifelse(is.na(id) | id == 0,
                       paste0("S", SERIAL, "_", PERNUM, "_", year), id)) %>%
    filter(!is.na(state))

  fwrite(data.cps2, file.path(DATA_GEN, "cps.csv"))
  cps_done <- TRUE
}

if (!cps_done) {
  # Fall back to pre-cleaned CSV (check Generated first, then Raw)
  cps_file <- file.path(DATA_GEN, "cps.csv")
  if (!file.exists(cps_file) || file.size(cps_file) == 0) {
    cps_file <- file.path(DATA_RAW, "cps.csv")
  }
  if (file.exists(cps_file) && file.size(cps_file) > 0) {
    cat("  Using pre-cleaned cps.csv\n")
    if (cps_file != file.path(DATA_GEN, "cps.csv")) {
      file.copy(cps_file, file.path(DATA_GEN, "cps.csv"), overwrite = TRUE)
    }
  } else {
    warning("No CPS data found. Skipping CPS cleaning.")
    cat("CPS cleaning skipped (no data).\n")
  }
}

cat("CPS cleaning complete.\n")
