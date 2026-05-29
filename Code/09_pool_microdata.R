# 09_pool_microdata.R — Pool CPS and PSID microdata

cat("Pooling microdata...\n")

# Read CPS
cps_file <- file.path(DATA_GEN, "cps.csv")
if (file.exists(cps_file) && file.info(cps_file)$size > 100) {
  cps.in <- fread(cps_file)

  # Select and standardize CPS variables
  cps_vars <- c("year", "id", "state", "r_77", "fipstate", "income", "wt", "age", "sex",
                "black", "self.employed", "ue", "emp", "health.shock", "fulltime",
                "rents", "poverty", "married", "divorced_sep", "pub_sub", "educ.yrs")
  available_cps <- intersect(cps_vars, names(cps.in))

  # Add owns if OWNERSHP exists
  if ("OWNERSHP" %in% names(cps.in)) {
    cps.in[, owns := ifelse(OWNERSHP == 22, 1, 0)]
  }
  if ("IND" %in% names(cps.in)) {
    cps.in[, industry := IND]
  }
  if ("METFIPS" %in% names(cps.in)) {
    cps.in[, fipmsa := METFIPS]
  }

  cps <- cps.in %>%
    dplyr::select(any_of(c(cps_vars, "owns", "industry", "fipmsa"))) %>%
    mutate(source = "CPS")
} else {
  # Try pre-built file
  pre_cps <- file.path(DATA_RAW, "cps.csv")
  if (file.exists(pre_cps)) {
    cps <- fread(pre_cps) %>% mutate(source = "CPS")
  } else {
    cps <- data.table()
    warning("No CPS data available for pooling")
  }
}

# Read PSID
psid_file <- file.path(DATA_GEN, "psid_final.csv")
if (file.exists(psid_file) && file.info(psid_file)$size > 100) {
  psid.in <- fread(psid_file)

  psid <- psid.in %>%
    mutate(rents = 1 - homeowner) %>%
    dplyr::select(any_of(c("year", "id", "State", "y1977", "income", "age", "black",
                            "rents", "owns", "ncars", "childcare_expense", "retirement",
                            "bankrupt", "move", "hh_size", "food_expense",
                            "food_insufficient", "insurance", "rent_expense"))) %>%
    rename(any_of(c(state = "State", r_77 = "y1977"))) %>%
    mutate(source = "PSID")
} else {
  psid <- data.table()
  warning("No PSID data available for pooling")
}

# Pool
if (nrow(cps) > 0 || nrow(psid) > 0) {
  pooled <- rbind(cps, psid, fill = TRUE)
  fwrite(pooled, file.path(DATA_GEN, "pooled_microdata.csv"))
  cat("Pooled microdata:", nrow(pooled), "observations\n")
} else {
  # Use pre-built pooled data if available
  pre_pooled <- file.path(DATA_RAW, "pooled_microdata.csv")
  if (file.exists(pre_pooled)) {
    file.copy(pre_pooled, file.path(DATA_GEN, "pooled_microdata.csv"), overwrite = TRUE)
    cat("Using pre-built pooled_microdata.csv\n")
  } else {
    fwrite(data.table(), file.path(DATA_GEN, "pooled_microdata.csv"))
    warning("No microdata available for pooling")
  }
}

cat("Microdata pooling complete.\n")
