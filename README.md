# Replication Package: Credit Cards and Retail Firms

## Overview

This replication package reproduces all figures, tables, the compiled paper, and presentation slides for:

**"Credit Cards and Retail Firms: Historical Evidence from the United States"**
by Joseph P. Hall (Georgia Institute of Technology)

## Requirements

### Software
- **R** (version 4.0 or later)
- **LaTeX** distribution with `latexmk` (e.g., TeX Live or MiKTeX)
- Required R packages are installed automatically by `Code/00_setup.R`

### Key R Packages
`data.table`, `dplyr`, `tidyr`, `ggplot2`, `fixest`, `haven`, `readxl`, `readr`, `stringr`, `sf`, `kableExtra`, `nleqslv`, `modelsummary`, `ipumsr`, `numDeriv`, `survey`

## How to Run

### One-click execution

```bash
bash run.sh
```

This will:
1. Install all required R packages
2. Run the full analysis pipeline (data cleaning + analysis)
3. Compile the paper PDF
4. Compile the slide deck PDF

### Step-by-step

```bash
# Install packages
Rscript Code/00_setup.R

# Run full pipeline
Rscript Code/run_all.R

# Compile paper
cd Paper && latexmk -pdf -interaction=nonstopmode main.tex
```

## Directory Structure

```
Replication/
├── run.sh                    # One-click execution script
├── Code/
│   ├── 00_setup.R            # Package installation and path configuration
│   ├── helpers/
│   │   ├── nlfunctions.R     # Nested logit demand model functions
│   │   └── scfcombo.R        # SCF multiple imputation variance estimation
│   │
│   ├── 01_clean_usury.R      # Usury rate data + choropleth maps
│   ├── 02_clean_scf.R        # Survey of Consumer Finances
│   ├── 03_clean_cbp.R        # County Business Patterns
│   ├── 04_clean_cps.R        # Current Population Survey (IPUMS)
│   ├── 05_clean_psid.R       # Panel Study of Income Dynamics
│   ├── 06_clean_bea.R        # Bureau of Economic Analysis
│   ├── 07_clean_compustat.R  # Compustat firm financials
│   ├── 08_clean_db.R         # Dun & Bradstreet establishment data
│   ├── 09_pool_microdata.R   # Pool CPS and PSID
│   │
│   ├── 10_summary_stats.R    # Summary statistics and time-series figures
│   ├── 11_did_bea_cbp.R      # DiD: BEA/CBP macro outcomes
│   ├── 12_did_compustat.R    # DiD: Compustat receivables/revenue
│   ├── 13_did_marquette.R    # DiD: D&B firm exit rates
│   ├── 14_did_pooled.R       # DiD: Household bankruptcy/self-employment
│   ├── 15_analyze_scranton.R # Credit card adoption (Scranton)
│   ├── 16_analyze_phonebooks.R # Yellow Pages card acceptance
│   ├── 17_analyze_receivables.R # QFR receivables analysis
│   ├── 18_analyze_occupations.R # Occupation composition trends
│   ├── 19_card_introduction.R   # Card introduction timing
│   ├── 20_permutation_test.R    # Permutation inference
│   ├── 21_nested_logit.R        # Nested logit demand model
│   ├── 22_macro_model.R         # General equilibrium model
│   └── run_all.R                # Master script (sources 00-22)
│
├── Data/
│   ├── Raw/              # Source data files (do not modify)
│   └── Generated/        # Intermediate and final datasets (created by pipeline)
│
├── Output/
│   ├── Figures/          # Generated figures (PDF/PNG)
│   └── Tables/           # Generated tables (TeX)
│
├── Paper/
│   ├── main.tex          # Paper source
│   ├── slides.tex        # Presentation slides (beamer)
│   ├── references.bib    # Bibliography
│   └── Static/           # Non-generated figures (scans, external)
│
└── README.md
```

## Data Sources

| Dataset | Source | Scripts |
|---------|--------|---------|
| Survey of Consumer Finances (SCF) | Federal Reserve | 02, 10 |
| County Business Patterns (CBP) | Census Bureau | 03, 10, 11 |
| Current Population Survey (CPS) | IPUMS-CPS | 04, 09, 14, 18 |
| Panel Study of Income Dynamics (PSID) | University of Michigan | 05, 09, 14 |
| Bureau of Economic Analysis (BEA) | BEA.gov | 06, 11 |
| Compustat | S&P Global via WRDS | 07, 12, 20 |
| Dun & Bradstreet | Historical establishment records | 08, 13, 15, 16 |
| Quarterly Financial Reports (QFR) | Census Bureau | 17 |
| Yellow Pages | Historical phonebooks | 16, 19 |
| Usury Rate Laws | Author compilation | 01 |
| FRED (Fed Funds, CPI) | Federal Reserve Bank of St. Louis | 01 |

## Output Mapping

### Figures
| Figure | File | Script |
|--------|------|--------|
| Spending levels | `spending_levels.pdf` | 10 |
| Store credit fraction | `store_fraction.pdf` | 10 |
| Receivables time series | `receivables.pdf` | 17 |
| Back-office vs retail | `backoffice_retail.pdf` | 18 |
| International comparison | `international.pdf` | 10 |
| Income bunching (1969) | `bunching.pdf` | 10 |
| Income bunching (1983) | `bunching_83.pdf` | 10 |
| Inflation time series | `inflation_ts.pdf` | 10 |
| Card growth | `card_growth_stacked.pdf` | 10 |
| Usury rate map (1977) | `rate_77.png` | 01 |
| Receivables/revenue DiD | `receivables_revenue.pdf` | 12 |
| BEA/CBP event study | `est.pdf` | 11 |
| Small establishments | `n1_19.pdf`, `frac_small.pdf` | 11 |
| Firm exit DiD | `difdif_exit_all.pdf` | 13 |
| Bankruptcy DiD | `bankruptcy.pdf` | 14 |
| Phonebooks by year | `books_by_year.pdf` | 19 |
| Card introduction | `intro_lb.png` | 19 |
| Interest rate histogram | `rates_all.png` | 10 |
| Permutation test | `permutation_test.pdf` | 20 |

### Tables
| Table | File | Script |
|-------|------|--------|
| SCF summary statistics | `summary_scf.tex` | 10 |
| CBP summary statistics | `summary_cbp.tex` | 10 |
| Card usage | `card_uses.tex` | 10 |
| Compustat regressions | `compustat_regressions.tex` | 12 |
| Firm adoption | `adoption_firms.tex` | 15 |
| Adoption effects | `adoption_effects.tex` | 16 |
| Industry adoption | `adoption_industries.tex` | 15 |
| Recession DiD | `recession.tex` | 13 |

## Notes

- Some D&B establishment-level data files are proprietary and may not be available. The pipeline includes fallback logic that uses pre-generated outputs when raw data is unavailable.
- All Stata code from the original project has been translated to R.
- The pipeline uses `fixest` for difference-in-differences estimation with two-way fixed effects.
