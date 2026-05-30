# Research Agenda — Archaeology of `/workspace/Input`

*Compiled May 2026. Source: 3-agent parallel sweep of `/workspace/Input` (old drafts,
slides, code, figures, tables) cross-checked against current data and code in
`/workspace/Replication/`. Purpose: catalog abandoned-but-promising analyses ("ghosts")
and novel threads worth reviving now that we have the bandwidth to pursue them.*

Status legend:
- **Runnable today** — data is already in `/workspace/Replication/`, just needs code.
- **Needs data** — analysis existed before; the input series must be relocated/rebuilt.
- **Novel** — would require genuinely new work; included deliberately ("ghosts are fine").

---

## A. Reduced-form mechanism evidence (data in hand)

### A1. Intensive-margin sales DiD — *Runnable today*
- **Question:** Conditional on survival, do treated retail firms *shrink* (lose sales) after
  Marquette, or is the action purely on the extensive (entry/exit) margin?
- **Old artifact:** `difdif_log_sls.pdf`
- **Approach:** Event-study DiD on `log_sls` / `sls_change` for the balanced panel of
  survivors. Mirror the exit specification.
- **Data:** `Data/Raw/DB/retail.csv` — has `log_sls`, `sls_change`, `sls_growth_dummy`,
  `sls_fall_dummy`, `treated`, `r_77`, `YEAR`, `DUNSNO`, `STATE`, SIC1–6.
- **Merit:** Cleanly separates extensive vs intensive margin; directly supports the
  entry-threshold story (B6 in the revision plan). Cheap.

### A2. Market structure / concentration (HHI) DiD — *Runnable today — TOP PICK*
- **Question:** Did entry induced by cheaper credit *deconcentrate* retail markets? This is
  the reduced-form counterpart of the structural markup/σ channel.
- **Old artifacts:** `Market Structure.pdf` (conceptual diagram), `simulation_small_share.png`
  (model comparative static: small-firm share rises with the rate cap).
- **Approach:** Compute Herfindahl (HHI) and small-firm sales share by SIC×state×year from
  firm-level `SLS`. Event-study DiD of HHI / small-share on treatment.
- **Data:** `retail.csv` (SLS, SIC, state, year, size flags `single`/`small`).
- **Merit:** Gives the markup/competition channel an **independent reduced-form leg**, which
  is exactly what R2 wanted (don't lean on the calibrated structural model for the
  competition result). Strongest single addition.

### A3. Size-heterogeneity gradient — *Runnable today*
- **Question:** Are small/single-location firms the ones that enter and grow, as the
  bank-dependence story predicts?
- **Old artifacts:** `dif_n20_49`, `growth_sizes`
- **Approach:** Split or interact treatment by employment-size bins / `single` / `small`.
- **Data:** `retail.csv` (`EMPLOYEESTHISSITE`, size codes, `single`, `small`).
- **Merit:** Tests the core mechanism (small firms depend on banks, large firms self-finance).

### A4. Industry heterogeneity — *Runnable today*
- **Question:** Does the effect concentrate in big-ticket / durable / clothing retail
  (where store credit mattered) and vanish in food (cash/grocery)?
- **Old artifact:** `food.pdf`
- **Approach:** Heterogeneity by SIC1 / 2-digit SIC; food (54xx) as a placebo.
- **Data:** `retail.csv` SIC fields.
- **Merit:** A credible placebo + dose-response across industries strengthens identification.

---

## B. Reduced-form, needs locating/rebuilding data

### B1. Price pass-through (state-CPI event study) — *Needs data — HIGHEST CEILING*
- **Question:** Did cheaper credit and tougher competition lower *retail prices* in treated
  states?
- **Old artifacts:** `price_level.pdf` (treated −1 log pt post-Marquette, flat pre-trend),
  `price_growth.pdf`
- **Approach:** State-level CPI event study, same DiD frame.
- **Data status:** A state-CPI series *existed* and produced a clean result. Likely BLS city
  CPI or Del Negro (1998) state-price reconstruction. Hunt in `Input/Data/inflation.dta`
  and `Input/Data/Local finances.zip`.
- **Merit:** A consumer-welfare/price result is the highest-value reduced-form outcome —
  speaks directly to the policy framing R2 requested. Worth the data hunt.

### B2. Tradable vs non-tradable decomposition (Mian–Sufi style) — *Needs data*
- **Question:** Is the employment/output response concentrated in non-tradable (local
  retail) vs tradable sectors?
- **Old artifacts:** `emp_tradable`, `emp_nontradable`, `gdp_tradable`, `gdp_nontradable`
- **Approach:** Classify SIC into tradable/non-tradable, run sectoral DiD on BEA employment/GDP.
- **Data status:** Needs an SIC→tradability crosswalk + sectoral BEA series.
- **Merit:** Connects to a well-known local-demand literature; moderate lift.

---

## C. Novel "ghost" threads (deliberate — may need new work)

### C3. Vertical-integration / double-marginalization theory — *Novel (theory)*
- **Idea:** `Market Structure.pdf` frames large retailers as *vertically integrated* across
  the credit market and the product market, while small retailers buy credit from banks
  (double marginalization). Marquette compresses the small-firm credit markup, shrinking the
  large-firm integration advantage.
- **Work:** Formalize as a vertical-relations model; test the implied markup-gap compression
  with A2's HHI/small-share evidence.
- **Merit:** Could become the paper's theoretical spine — reframes the structural model as a
  test of a specific IO mechanism rather than a generic nested logit.

### C1. Geographic distance-decay of competition — *Novel (own paper potential)*
- **Idea:** Old `estimate_t.R` geocoded D&B street addresses and found ~4% less sales decline
  per 1% more distance to the nearest entrant.
- **Work:** Re-geocode `STREETADDRESS`/`CITY`/`ZIPCODE` in retail.csv; build distance-to-entrant.
- **Merit:** A spatial-competition result; potentially a standalone paper.

### C2. Bank-side risk absorption — *Novel — BoC-relevant*
- **Question:** Did banks in treated states absorb more consumer-credit losses in the early
  1980s as they expanded card lending?
- **Data status:** `Input/Data/Call reports.zip` and `Input/Data/FFIEC.zip` are present.
- **Hard part:** Mapping borrower location to bank location (banks lend across state lines
  post-Marquette by design).
- **Merit:** Financial-stability angle — directly relevant to the Bank of Canada audience.

### C4. Own-credit advertising retreat — *Novel*
- **Idea:** ~120k Yellow Pages "own credit" ads were parsed; DiD on the dismantling of
  in-house store-credit promotion after bank cards arrive.
- **Merit:** Novel mechanism evidence using the YP corpus we already process.

### C5. Back-office labor displacement DiD — *Novel (prior null)*
- **Idea:** Old `did_occupations.R` (IPUMS) tested whether bank cards displaced store credit
  back-office clerical labor; prior runs were null.
- **Merit:** Even a tight null is an informative bound; low cost if IPUMS extract survives.

### C6. Consumer-side balance sheets — *Novel*
- **Idea:** `networth_negative`, `bankruptcy` outcomes from PSID/SCF.
- **Data status:** `Input/Data/PSID.zip`, `Input/Data/SCF.zip` present.
- **Merit:** Household-welfare counterpart to the firm-side story.

### C7. Entry composition (not just net entry) — *Runnable-ish*
- **Question:** *What kind* of firms entered — size, industry, single vs chain?
- **Data:** retail.csv `entry` flag + firm characteristics.
- **Merit:** Refines the headline entry result; cheap.

---

## D. Presentation / framing assets (for BoC talk + R2 policy ask)
- Afterpay / Square / BNPL images, UK regulation slides, Chase / Walmart / imprinter
  historicals, `intro_exact` + `graphics1–4` mockups.
- **Use:** Motivate present-day relevance (interchange / BNPL regulation), satisfying R2's
  policy-relevance request and sharpening the Bank of Canada framing.

---

## E. Correctly dead (do not revive)
- Modern Experian / RateWatch ML heterogeneity / pass-through / credit-limit dynamics —
  belongs to the separate `jphall-loans` project; proprietary data; era-mismatch.
- 1958 Fresno BankAmericard pilot — entry was endogenous; not a clean experiment.
- Auto-loan analyses — tangential to the retail-firm question.

---

## Suggested sequencing
1. **A2 (HHI/market structure)** + **A1 (intensive margin)** + **A4 (industry placebo)** —
   all runnable today from retail.csv; together they give the competition channel an
   independent reduced-form leg.
2. **B1 (price pass-through)** — locate the state-CPI series; highest-value outcome.
3. **C3 (vertical-integration theory)** to tie A2's evidence to a clean mechanism.
4. Opportunistic: C2 (bank-side, BoC-relevant), C7 (entry composition), D (framing assets).

*Resume agent IDs (if deeper digging needed): drafts `a2dcf76fb6bc4d863`,
code `a1920eab1adbbfd09`, figures `acd6c73bb7e5ebabb`.*
