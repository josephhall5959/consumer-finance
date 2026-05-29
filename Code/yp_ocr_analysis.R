# yp_ocr_analysis.R — Analyze raw OCR data from Yellow Pages phonebooks
# 110 CSV files from Hi-Tech Inc. OCR processing of Library of Congress scans

cat("Analyzing Yellow Pages OCR data...\n")

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

OCR_DIR <- "/workspace/Input/_extracted/Yellow_Pages_OCR/Merged Files"
OUTPUT <- "Output/Figures/yp_review"

# ============================================================
# Pass 1: Lightweight scan of all 110 files
# ============================================================
cat("  Pass 1: Scanning all phonebook files...\n")

files <- list.files(OCR_DIR, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
cat(sprintf("  Found %d files\n", length(files)))

book_summaries <- rbindlist(lapply(seq_along(files), function(i) {
  f <- files[i]
  tryCatch({
    # Read full file but only needed columns
    d <- fread(f, select = c("index", "set", "state", "city_county", "year",
                              "page", "is_title", "title_corrected",
                              "has.number", "text", "body.height"),
               showProgress = FALSE)

    # Extract book ID from filename
    book_id <- gsub(".*book_([0-9]+).*", "\\1", basename(f))

    # Count titles (Yellow Pages category headers)
    titles <- d[is_title == 1 & !is.na(title_corrected)]
    unique_categories <- unique(titles$title_corrected)
    # Clean up category names
    unique_categories <- unique_categories[nchar(unique_categories) > 1]
    unique_categories <- unique(gsub("\\s*\\(Cont'd\\)$|\\s*-\\s*$", "",
                                     unique_categories, ignore.case = TRUE))

    # Count listings (rows that look like business entries: have phone numbers)
    n_listings <- sum(d$has.number == 1, na.rm = TRUE)

    # Pages
    pages <- unique(d$page)

    data.table(
      book_id = book_id,
      file = basename(f),
      state = d$state[1],
      city = d$city_county[1],
      year = d$year[1],
      n_rows = nrow(d),
      n_pages = length(pages),
      n_listings = n_listings,
      n_categories = length(unique_categories),
      n_sets = length(unique(d$set)),
      categories_sample = paste(head(sort(unique_categories), 10), collapse = "; ")
    )
  }, error = function(e) {
    data.table(book_id = basename(f), file = basename(f),
               state = NA, city = NA, year = NA,
               n_rows = 0, n_pages = 0, n_listings = 0,
               n_categories = 0, n_sets = 0, categories_sample = "ERROR")
  })
}))

cat(sprintf("  Scanned %d phonebooks\n", nrow(book_summaries)))
cat(sprintf("  Total OCR words: %s\n", format(sum(book_summaries$n_rows), big.mark = ",")))
cat(sprintf("  Total pages: %s\n", format(sum(book_summaries$n_pages), big.mark = ",")))
cat(sprintf("  Total listings (with phone numbers): %s\n",
            format(sum(book_summaries$n_listings), big.mark = ",")))
cat(sprintf("  States: %s\n", paste(unique(book_summaries$state), collapse = ", ")))
cat(sprintf("  Year range: %s-%s\n", min(book_summaries$year, na.rm = TRUE),
            max(book_summaries$year, na.rm = TRUE)))

# ============================================================
# Pass 2: Extract ALL category headers across all phonebooks
# ============================================================
cat("  Pass 2: Extracting category headers...\n")

all_categories <- rbindlist(lapply(files, function(f) {
  tryCatch({
    d <- fread(f, select = c("state", "city_county", "year", "page",
                              "is_title", "title_corrected"),
               showProgress = FALSE)
    titles <- d[is_title == 1 & !is.na(title_corrected) & nchar(title_corrected) > 1]
    if (nrow(titles) == 0) return(NULL)
    titles[, .(state = state[1], city = city_county[1], year = year[1],
               category = gsub("\\s*\\(Cont'd\\)$|\\s*-\\s*$|\\s*\\(Cont.d\\)$", "",
                               title_corrected, ignore.case = TRUE),
               page = page)]
  }, error = function(e) NULL)
}))

# Clean and count
all_categories[, category_clean := str_to_title(trimws(category))]
cat_counts <- all_categories[, .N, by = category_clean][order(-N)]
cat(sprintf("  Found %d unique category names across all phonebooks\n", nrow(cat_counts)))

# ============================================================
# Figure 7: OCR Data Summary — Phonebooks by state and year
# ============================================================
cat("  Creating OCR summary figures...\n")

book_summaries[, year_fct := factor(year)]
p7 <- ggplot(book_summaries[!is.na(year)], aes(x = year_fct, fill = state)) +
  geom_bar() +
  theme_classic(base_size = 13) +
  labs(x = "Year", y = "Number of phonebooks (OCR-processed)",
       fill = "State") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT, "7_ocr_books_by_year.pdf"), p7, width = 9, height = 5)
cat("    Created 7_ocr_books_by_year.pdf\n")

# ============================================================
# Figure 8: Distribution of listings per phonebook
# ============================================================
p8 <- ggplot(book_summaries, aes(x = n_listings)) +
  geom_histogram(bins = 25, fill = "steelblue", alpha = 0.8, color = "white") +
  theme_classic(base_size = 13) +
  labs(x = "Business listings per phonebook (entries with phone numbers)",
       y = "Number of phonebooks") +
  geom_vline(xintercept = median(book_summaries$n_listings),
             linetype = "dashed", color = "red") +
  annotate("text", x = median(book_summaries$n_listings), y = Inf,
           label = paste0("Median = ", format(median(book_summaries$n_listings), big.mark = ",")),
           vjust = 2, hjust = -0.1, color = "red", size = 4)
ggsave(file.path(OUTPUT, "8_listings_per_book.pdf"), p8, width = 8, height = 5)
cat("    Created 8_listings_per_book.pdf\n")

# ============================================================
# Figure 9: Top 30 most common business categories
# ============================================================
top_cats <- head(cat_counts, 30)
p9 <- ggplot(top_cats, aes(x = reorder(category_clean, N), y = N)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Appearances across all phonebooks") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(file.path(OUTPUT, "9_top_categories.pdf"), p9, width = 8, height = 7)
cat("    Created 9_top_categories.pdf\n")

# ============================================================
# Figure 10: Pages per phonebook over time
# ============================================================
book_summaries[, year_num := as.numeric(year)]
p10 <- ggplot(book_summaries[!is.na(year_num)], aes(x = year_num, y = n_pages)) +
  geom_point(aes(color = state), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  theme_classic(base_size = 13) +
  labs(x = "Year", y = "Pages per phonebook (OCR-processed)",
       color = "State") +
  scale_color_brewer(palette = "Set2")
ggsave(file.path(OUTPUT, "10_pages_over_time.pdf"), p10, width = 8, height = 5)
cat("    Created 10_pages_over_time.pdf\n")

# ============================================================
# Figure 11: Credit-card-related categories
# ============================================================
credit_cats <- all_categories[grepl("credit|charge|bank.?americ|master.?charge|visa",
                                     category, ignore.case = TRUE)]
if (nrow(credit_cats) > 0) {
  credit_summary <- credit_cats[, .(n_entries = .N), by = .(city, state, year)]
  credit_summary <- credit_summary[order(year)]

  p11 <- ggplot(credit_summary, aes(x = year, y = n_entries)) +
    geom_col(aes(fill = state), alpha = 0.8) +
    theme_classic(base_size = 13) +
    labs(x = "Year", y = "Credit-related category entries in Yellow Pages",
         fill = "State") +
    scale_fill_brewer(palette = "Set2")
  ggsave(file.path(OUTPUT, "11_credit_categories_over_time.pdf"), p11, width = 8, height = 5)
  cat("    Created 11_credit_categories_over_time.pdf\n")

  cat("\n  Credit-related categories found:\n")
  print(credit_cats[, .N, by = category_clean][order(-N)])
} else {
  cat("    No credit-related categories found in title headers\n")
}

# ============================================================
# Summary statistics table
# ============================================================
ocr_summary <- data.frame(
  Metric = c(
    "Phonebooks processed (OCR)",
    "Total OCR words extracted",
    "Total pages scanned",
    "Business listings (with phone numbers)",
    "Unique business categories",
    "States covered",
    "Cities covered",
    "Year range",
    "Median listings per phonebook",
    "Median pages per phonebook",
    "Median categories per phonebook"
  ),
  Value = c(
    nrow(book_summaries),
    format(sum(book_summaries$n_rows), big.mark = ","),
    format(sum(book_summaries$n_pages), big.mark = ","),
    format(sum(book_summaries$n_listings), big.mark = ","),
    format(nrow(cat_counts), big.mark = ","),
    length(unique(book_summaries$state)),
    length(unique(book_summaries$city)),
    paste0(min(book_summaries$year, na.rm = TRUE), "-",
           max(book_summaries$year, na.rm = TRUE)),
    format(median(book_summaries$n_listings), big.mark = ","),
    median(book_summaries$n_pages),
    median(book_summaries$n_categories)
  )
)

# Save as LaTeX
library(kableExtra)
tex_out <- ocr_summary %>%
  kbl(format = "latex", booktabs = TRUE, col.names = c("", "")) %>%
  kable_classic() %>%
  as.character()
tex_lines <- strsplit(tex_out, "\n")[[1]]
tex_lines <- tex_lines[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption",
                               tex_lines)]
writeLines(tex_lines, file.path(OUTPUT, "ocr_summary.tex"))
cat("  Created ocr_summary.tex\n")

# Save book-level summaries for reference
fwrite(book_summaries, file.path(OUTPUT, "book_summaries.csv"))
cat("  Created book_summaries.csv\n")

cat("\nOCR analysis complete.\n")
