# 18_analyze_occupations.R — Analyze occupational composition changes
# Based on original "Analyze Occupations.R"
# Uses IPUMS USA decennial Census microdata (usa_00008)

cat("Analyzing occupations...\n")

library(ipumsr)

ddi_file <- file.path(DATA_RAW, "CPS", "usa_00008.xml")
dat_file <- file.path(DATA_RAW, "CPS", "usa_00008.dat")

if (!file.exists(ddi_file) || !file.exists(dat_file) || file.size(dat_file) < 1000) {
  cat("  usa_00008 data not found or empty. Skipping occupation analysis.\n")
} else {
  # Read fixed-width file directly using column positions from DDI
  # Columns: YEAR(1-4), PERWT(82-91), EMPSTAT(92), OCC1950(99-101), IND1950(102-104), INCWAGE(105-110)
  col_positions <- fwf_positions(
    start = c(1, 82, 92, 99, 102, 105),
    end   = c(4, 91, 92, 101, 104, 110),
    col_names = c("YEAR", "PERWT", "EMPSTAT", "OCC1950", "IND1950", "INCWAGE")
  )
  cat("  Reading fixed-width data in chunks to conserve memory...\n")
  chunk_size <- 5e6
  n_lines <- as.integer(system(paste("wc -l <", dat_file), intern = TRUE))
  cat(sprintf("  Total lines: %d, reading in chunks of %d\n", n_lines, chunk_size))

  ts_list <- list()
  offset <- 0
  chunk_i <- 1
  while (offset < n_lines) {
    chunk <- read_fwf(dat_file, col_positions,
                       col_types = "diiiid",
                       skip = offset, n_max = chunk_size,
                       progress = FALSE) %>%
      as.data.table()
    chunk <- chunk[YEAR >= 1940]
    chunk[, employed := fifelse(EMPSTAT == 1L, 1L, 0L)]
    chunk[, wholesale := fifelse(IND1950 >= 606L & IND1950 <= 627L, 1L, 0L)]
    chunk[, retail := fifelse(IND1950 >= 636L & IND1950 <= 699L, 1L, 0L)]
    chunk[, back_office := fifelse(OCC1950 %in% c(310L, 321L), 1L, 0L)]
    chunk[, wage := fifelse(INCWAGE == 999999, NA_real_, as.numeric(INCWAGE))]

    ts_chunk <- chunk[, .(
      wages_retail = sum(wage * employed * retail * PERWT / 100, na.rm = TRUE),
      wages_wholesale = sum(wage * employed * wholesale * PERWT / 100, na.rm = TRUE),
      wages_retail_backoffice = sum(wage * employed * retail * back_office * PERWT / 100, na.rm = TRUE),
      wages_wholesale_backoffice = sum(wage * employed * wholesale * back_office * PERWT / 100, na.rm = TRUE)
    ), by = YEAR]
    ts_list[[chunk_i]] <- ts_chunk
    rm(chunk); gc(verbose = FALSE)
    offset <- offset + chunk_size
    chunk_i <- chunk_i + 1
  }
  cat(sprintf("  Processed %d chunks\n", chunk_i - 1))

  ts <- rbindlist(ts_list)[, .(
    wages_retail = sum(wages_retail),
    wages_wholesale = sum(wages_wholesale),
    wages_retail_backoffice = sum(wages_retail_backoffice),
    wages_wholesale_backoffice = sum(wages_wholesale_backoffice)
  ), by = YEAR]
  rm(ts_list); gc(verbose = FALSE)
  ts[, backoffice_in_retail := wages_retail_backoffice / wages_retail]
  ts[, backoffice_in_wholesale := wages_wholesale_backoffice / wages_wholesale]
  setorder(ts, YEAR)

  cat("  Backoffice wage shares:\n")
  print(ts[, .(YEAR, backoffice_in_retail, backoffice_in_wholesale)])

  # Figure: Bookkeepers and bill collectors share of retail wage bill
  p <- ggplot(ts, aes(x = YEAR, y = backoffice_in_retail)) +
    geom_line() +
    labs(y = "Share of Retail Wage Bill",
         title = "Bookkeepers and Bill Collectors") +
    theme_classic()
  ggsave(file.path(OUTPUT_FIG, "backoffice_retail.pdf"), p, width = 9, height = 6)
  cat("  Created backoffice_retail.pdf\n")

  # Report drop statistics
  if ("1970" %in% as.character(ts$YEAR)) {
    val_1970 <- ts[YEAR == 1970, backoffice_in_retail]
    val_last <- ts[.N, backoffice_in_retail]
    val_w_1970 <- ts[YEAR == 1970, backoffice_in_wholesale]
    val_w_last <- ts[.N, backoffice_in_wholesale]
    cat(sprintf("  Retail drop since 1970: %.0f%%\n", 100 * (1 - val_last / val_1970)))
    cat(sprintf("  Wholesale drop since 1970: %.0f%%\n", 100 * (1 - val_w_last / val_w_1970)))
  }

  rm(data.in, ts); gc()
}

cat("Occupation analysis complete.\n")
