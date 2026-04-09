################################################################################
##
##  Synthetic ICU Extubation Cohort — Data Generator
##
##  PURPOSE
##  -------
##  Generates a synthetic patient-level dataset structurally analogous to the
##  MIMIC-IV last_extubations cohort used in MIMIC_analysis_clean.R, for use
##  in public-facing visualizations (e.g., Tableau Public dashboards).
##
##  NO MIMIC DATA ARE USED HERE. All distributional parameters are sourced from
##  published literature (citations inline). The simulation preserves:
##    - Realistic marginal distributions for each variable
##    - Clinically plausible correlations among predictors
##    - Outcome rates consistent with published ICU extubation literature
##    - The caregiver-level structure (FE rate, volume) central to the analysis
##
##  This approach is methodologically appropriate for portfolio/demonstration
##  purposes and is defensible in interviews: synthetic generation preserving
##  distributional structure is a standard technique when data use agreements
##  preclude sharing patient-level records.
##
##  PARAMETER SOURCES
##  -----------------
##  Age:
##    MIMIC-III published demographics: median 65.8, IQR 52.8–77.8
##    (Johnson et al., Scientific Data 2016)
##
##  Failed extubation rate (overall):
##    Published range 10–20% across ICU types; 15.4% in Fernandez et al. 2024
##    (Respiratory Care); 9–17% in surgical ICU (PubMed 34756824)
##
##  12-month survival:
##    Hospital mortality ~11.5% in MIMIC-III broadly; extubation failure group
##    42% in-hospital mortality vs 14% in success group (Fernandez et al. 2024)
##    12-month survival modeled as somewhat lower than in-hospital survival
##
##  Charlson comorbidity index:
##    Median ~5 in MIMIC-IV ICU sepsis cohort (MIMIC-Sepsis, arXiv 2025)
##
##  SOFA score:
##    Moderate-to-severe illness in ventilated patients; mean ~7–9 in published
##    mechanical ventilation cohorts
##
##  Ventilation hours:
##    Median ~80–120h in general ICU vent cohorts; right-skewed (lognormal)
##    Filtered to (1, 1000) matching SQL pipeline
##
##  Norepinephrine equivalent dose:
##    Vasopressor use ~17–35% of ICU patients; filtered to <1 mcg/kg/min
##    matching SQL pipeline
##
##  Caregiver FE rate:
##    Simulated from Beta distribution; inter-caregiver variation documented
##    in nursing/provider outcome literature
##
##  ICD chapter prevalence:
##    ICU admission mix approximated from MIMIC-III/IV publications:
##    Circulatory ~30%, Respiratory ~20%, Digestive ~10%, Neurologic ~8%, etc.
##
##  OUTPUTS
##  -------
##  synthetic_last_extubations.csv  — patient-level data (N = 800)
##    Matches column structure expected by MIMIC_analysis_clean.R Section 2+
##  synthetic_caregiver_rate.csv    — caregiver-level summary (N ~ 80)
##    Pre-aggregated for Tableau caregiver views
##
################################################################################


## ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(MASS)       # mvrnorm for correlated predictors (load before dplyr to avoid select() conflict)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(stringr)
library(forcats)

set.seed(237)

N_PATIENTS   <- 800   # patient rows
N_CAREGIVERS <- 80    # unique caregivers (avg ~10 patients each, range 10–60)


## ── 1. CAREGIVER POOL ────────────────────────────────────────────────────────
##
##  Generate caregivers first, assign each a true FE rate drawn from a
##  Beta distribution. Volume (n_extubations) is drawn independently,
##  with mild positive skew (some high-volume nurses/residents).
##
##  Beta(2, 12): mean ~0.14, SD ~0.09 → plausible inter-caregiver spread
##  around the published 10–20% overall FE rate.

caregiver_pool <- tibble(
  caregiver_id     = sprintf("CG%04d", seq_len(N_CAREGIVERS)),
  true_fe_rate     = rbeta(N_CAREGIVERS, shape1 = 2, shape2 = 12),
  caregiver_n      = round(rlnorm(N_CAREGIVERS, meanlog = 3.0, sdlog = 0.7))
) %>%
  # Enforce caregiver_n > 10 (mirrors SQL filter caregiver_n > 10)
  mutate(caregiver_n = pmax(caregiver_n, 11)) %>%
  # Cap at 200 to avoid unrealistically high-volume single caregivers
  mutate(caregiver_n = pmin(caregiver_n, 200))


## ── 2. PATIENT-LEVEL PREDICTORS ──────────────────────────────────────────────
##
##  Continuous predictors are drawn from a multivariate normal then
##  back-transformed to their natural scales. Correlation structure
##  encodes known clinical relationships:
##    - charlson ↔ sofa: moderate positive (sicker patients have more comorbidities)
##    - charlson ↔ vent_hours: weak positive (more comorbid → longer ventilation)
##    - sofa ↔ norepinephrine: moderate positive (organ failure → vasopressors)
##    - sofa ↔ vent_hours: moderate positive (sicker → longer vent)
##    - age ↔ charlson: moderate positive (older → more comorbidities)
##    - age ↔ sofa: weak positive

Sigma <- matrix(c(
#         age    charlson  sofa   norepi  vent_log
          1.00,  0.40,     0.20,  0.10,   0.15,   # age
          0.40,  1.00,     0.35,  0.15,   0.25,   # charlson
          0.20,  0.35,     1.00,  0.45,   0.40,   # sofa
          0.10,  0.15,     0.45,  1.00,   0.20,   # norepinephrine
          0.15,  0.25,     0.40,  0.20,   1.00    # log(vent_hours)
), nrow = 5)

Z <- mvrnorm(
  n   = N_PATIENTS,
  mu  = rep(0, 5),
  Sigma = Sigma
)

patients_raw <- tibble(
  z_age      = Z[, 1],
  z_charlson = Z[, 2],
  z_sofa     = Z[, 3],
  z_norepi   = Z[, 4],
  z_vent_log = Z[, 5]
) %>%
  mutate(
    # Age: median 65, SD ~14, clipped to [18, 95]
    # Source: MIMIC-III median 65.8, IQR 52.8–77.8
    anchor_age  = pmin(pmax(round(65 + 14 * z_age), 18), 95),

    # Charlson comorbidity index: median ~5, right-skewed, clipped to [0, 15]
    # Source: MIMIC-Sepsis arXiv 2025 median CCI = 5
    charlson    = pmin(pmax(round(5 + 2.5 * z_charlson), 0), 15),

    # SOFA: mean ~8 in ventilated patients, SD ~3, clipped to [0, 20]
    sofa        = pmin(pmax(round(8 + 3 * z_sofa), 0), 20),

    # Norepinephrine equivalent: ~30% receive vasopressors; bounded (0, 1)
    # Zero-inflated: ~70% get 0, remainder from lognormal
    norepi_raw  = pmax(0.05 + 0.3 * z_norepi, 0.001),
    norepinephrine = ifelse(runif(N_PATIENTS) < 0.70, 0,
                            pmin(norepi_raw, 0.99)),

    # Vent hours: lognormal, median ~96h (~4 days), filtered (1, 1000)
    # Source: published ICU vent cohort medians ~80–120h
    vent_hours  = pmin(pmax(exp(4.6 + 0.8 * z_vent_log), 1), 999)
  ) %>%
  select(anchor_age, charlson, sofa, norepinephrine, vent_hours)


## ── 3. CATEGORICAL PREDICTORS ────────────────────────────────────────────────
##
##  gender: ~45% female in MIMIC ICU populations
##  intubation_type: surgical ~35%, medical-respiratory ~35%, medical-non-resp ~30%
##    (approximate from MIMIC-IV ICU composition and ICD chapter mix)
##  tube_event_source: ~60% explicit, ~40% inferred (reflects pipeline design)

patients_cat <- tibble(
  gender = sample(
    c("F", "M"),
    size    = N_PATIENTS,
    replace = TRUE,
    prob    = c(0.45, 0.55)
  ),
  intubation_type = sample(
    c("surgical", "medical-respiratory", "medical-non-respiratory"),
    size    = N_PATIENTS,
    replace = TRUE,
    prob    = c(0.35, 0.35, 0.30)
  ) %>% factor(),
  tube_event_source = sample(
    c("explicit", "inferred"),
    size    = N_PATIENTS,
    replace = TRUE,
    prob    = c(0.60, 0.40)
  ) %>% factor()
)


## ── 4. ICD CHAPTER ASSIGNMENT ────────────────────────────────────────────────
##
##  Prevalence approximated from MIMIC-III/IV ICU admission ICD distributions.
##  Constrained so intubation_type is correlated with ICD chapter:
##    - surgical → higher Circulatory, lower Respiratory
##    - medical-respiratory → higher Respiratory, Infectious
##    - medical-non-respiratory → higher Neurologic, Digestive, Circulatory

icd_chapters <- c(
  "Circulatory", "Respiratory", "Digestive", "Neurologic",
  "Infectious", "Genitourinary", "Musculoskeletal", "Pregnancy",
  "Symptoms", "SpecialCodes_e.g.COVID", "Endocrine", "Other"
)

assign_icd <- function(intub_type) {
  probs <- case_when(
    intub_type == "surgical" ~ list(c(
      0.38, 0.12, 0.14, 0.08, 0.06, 0.06, 0.05, 0.02, 0.03, 0.02, 0.02, 0.02
    )),
    intub_type == "medical-respiratory" ~ list(c(
      0.18, 0.32, 0.08, 0.07, 0.14, 0.05, 0.03, 0.02, 0.04, 0.04, 0.01, 0.02
    )),
    intub_type == "medical-non-respiratory" ~ list(c(
      0.28, 0.10, 0.14, 0.16, 0.09, 0.07, 0.04, 0.03, 0.04, 0.02, 0.02, 0.01
    ))
  )[[1]]
  sample(icd_chapters, size = 1, prob = probs)
}

primary_diagnosis <- map_chr(patients_cat$intubation_type, assign_icd) %>% factor()


## ── 5. CAREGIVER ASSIGNMENT ───────────────────────────────────────────────────
##
##  Assign patients to caregivers with roughly lognormal patient loads.
##  caregiver_fe_rate is the caregiver's true rate (used to generate outcomes).
##  caregiver_n is the caregiver's total volume (from pool, not patient count here).

caregiver_assignment <- sample(
  seq_len(N_CAREGIVERS),
  size    = N_PATIENTS,
  replace = TRUE,
  prob    = caregiver_pool$caregiver_n / sum(caregiver_pool$caregiver_n)
)

caregiver_fe_rate <- caregiver_pool$true_fe_rate[caregiver_assignment]
caregiver_n       <- caregiver_pool$caregiver_n[caregiver_assignment]
caregiver_id      <- caregiver_pool$caregiver_id[caregiver_assignment]


## ── 6. OUTCOME GENERATION ────────────────────────────────────────────────────
##
##  Outcomes are generated from a logistic model whose structure mirrors
##  logreg_icd in MIMIC_analysis_clean.R, with coefficients chosen to
##  reproduce published effect directions and plausible outcome rates.
##
##  Coefficients are rounded/stylized — they do not reproduce the actual
##  fitted coefficients from the MIMIC analysis.
##
##  TARGET RATES (from literature):
##    failed_extubation: ~15% overall
##    survival_12mo: ~80% (hospital mortality ~11%, additional 9% die within 12mo)
##
##  Key effects encoded:
##    - Higher caregiver_fe_rate → lower survival (quadratic)
##    - Higher charlson → lower survival
##    - Higher sofa → lower survival
##    - Higher norepinephrine → lower survival (strong)
##    - Longer vent_hours → modestly lower survival
##    - surgical intubation → better survival than medical
##    - Neurologic, Infectious, Digestive diagnoses → lower survival
##    - Circulatory, Respiratory diagnoses → relatively better survival

# Centered predictors (mirrors Section 2 of analysis script)
mean_charlson    <- mean(patients_raw$charlson)
mean_sofa        <- mean(patients_raw$sofa)
mean_norepi      <- mean(patients_raw$norepinephrine)
mean_vent        <- mean(patients_raw$vent_hours)

charlson_c    <- patients_raw$charlson     - mean_charlson
sofa_c        <- patients_raw$sofa        - mean_sofa
norepi_c      <- patients_raw$norepinephrine - mean_norepi
vent_c        <- patients_raw$vent_hours   - mean_vent

# ICD chapter risk adjustments (log-odds)
icd_effect <- case_when(
  primary_diagnosis == "Neurologic"             ~  0.40,
  primary_diagnosis == "Infectious"             ~  0.25,
  primary_diagnosis == "Digestive"              ~  0.20,
  primary_diagnosis == "SpecialCodes_e.g.COVID" ~  0.30,
  primary_diagnosis == "Circulatory"            ~ -0.20,
  primary_diagnosis == "Respiratory"            ~ -0.15,
  primary_diagnosis == "Genitourinary"          ~ -0.25,
  primary_diagnosis == "Pregnancy"              ~ -0.50,
  primary_diagnosis == "Musculoskeletal"        ~ -0.10,
  TRUE                                          ~  0.00
)

intub_effect <- case_when(
  patients_cat$intubation_type == "surgical"               ~ -0.30,
  patients_cat$intubation_type == "medical-respiratory"    ~  0.10,
  patients_cat$intubation_type == "medical-non-respiratory" ~  0.20
)

# Linear predictor for survival_12mo
# Intercept chosen so baseline survival ~80%
lp_survival <- 1.70 +
  intub_effect +
  icd_effect +
  (-0.18) * charlson_c +
  (-0.20) * sofa_c +
  (-1.80) * norepi_c +
  (-0.004) * vent_c +
  (-3.0)  * caregiver_fe_rate +
  ( 5.0)  * caregiver_fe_rate^2 +   # quadratic: very high FE rate caregivers
  rnorm(N_PATIENTS, sd = 0.5)        # residual patient-level noise

p_survival    <- plogis(lp_survival)
survival_12mo <- rbinom(N_PATIENTS, 1, p_survival) == 1

# Failed extubations: count (0–3), driven by caregiver FE rate + patient factors
# Baseline ~15% have at least one failed extubation
lp_fe <- -2.6 +
  (3.5) * caregiver_fe_rate +
  (0.08) * sofa_c +
  (0.05) * charlson_c +
  (0.80) * (patients_cat$intubation_type == "medical-non-respiratory") +
  rnorm(N_PATIENTS, sd = 0.4)

p_fe              <- plogis(lp_fe)
failed_extubation_any <- rbinom(N_PATIENTS, 1, p_fe)
failed_extubations    <- failed_extubation_any *
  sample(1:3, N_PATIENTS, replace = TRUE, prob = c(0.75, 0.20, 0.05))

cat("Simulated failed extubation rate:", mean(failed_extubations > 0), "\n")
cat("Simulated 12-month survival rate:", mean(survival_12mo), "\n")


## ── 7. ASSEMBLE PATIENT DATASET ──────────────────────────────────────────────

synthetic_patients <- bind_cols(
  tibble(
    subject_id        = sprintf("SYN%05d", seq_len(N_PATIENTS)),
    caregiver_id      = caregiver_id,
    caregiver_n       = caregiver_n,
    caregiver_fe_rate = caregiver_fe_rate
  ),
  patients_raw,
  patients_cat,
  tibble(
    primary_diagnosis  = primary_diagnosis,
    failed_extubations = failed_extubations,
    survival_12mo      = survival_12mo,
    hospital_expire_flag = as.integer(!survival_12mo & runif(N_PATIENTS) < 0.6)
  )
) %>%
  # Centered columns (mirrors cleandata in analysis script)
  mutate(
    vent_hours_centered     = vent_hours     - mean(vent_hours),
    charlson_centered       = charlson       - mean(charlson),
    norepinephrine_centered = norepinephrine - mean(norepinephrine),
    sofa_centered           = sofa           - mean(sofa)
  )

glimpse(synthetic_patients)
summary(synthetic_patients)


## ── 8. CAREGIVER-LEVEL SUMMARY ───────────────────────────────────────────────
##
##  Pre-aggregated for Tableau caregiver views.
##  caregiver_fe_rate here is computed from patient records (as in the analysis),
##  which will differ slightly from true_fe_rate due to sampling noise.

synthetic_caregiver_rate <- synthetic_patients %>%
  group_by(caregiver_id) %>%
  summarise(
    caregiver_n            = max(caregiver_n),
    caregiver_fe_rate      = max(caregiver_fe_rate),    # assigned rate
    observed_fe_rate       = mean(failed_extubations > 0),  # sample realization
    hospital_expire_rate   = mean(hospital_expire_flag),
    survival_12mo_rate     = mean(survival_12mo),
    n_patients_in_sample   = n(),
    .groups = "drop"
  )

cat("Caregiver summary:\n")
summary(synthetic_caregiver_rate)


## ── 9. VALIDATION CHECKS ─────────────────────────────────────────────────────
##
##  Verify simulated rates are plausible relative to published benchmarks.

cat("\n── VALIDATION ──────────────────────────────────────────────────\n")
cat("Overall FE rate (target ~10–20%):", round(mean(synthetic_patients$failed_extubations > 0), 3), "\n")
cat("Overall 12-mo survival (target ~75–85%):", round(mean(synthetic_patients$survival_12mo), 3), "\n")
cat("Median age (target ~63–67):", median(synthetic_patients$anchor_age), "\n")
cat("Median Charlson (target ~4–6):", median(synthetic_patients$charlson), "\n")
cat("Median SOFA (target ~7–9):", median(synthetic_patients$sofa), "\n")
cat("Median vent_hours (target ~80–120h):", round(median(synthetic_patients$vent_hours), 1), "\n")
cat("Vasopressor use (target ~25–35%):", round(mean(synthetic_patients$norepinephrine > 0), 3), "\n")
cat("Caregiver FE rate range:", round(range(synthetic_caregiver_rate$caregiver_fe_rate), 3), "\n")
cat("─────────────────────────────────────────────────────────────────\n\n")

# Intubation type breakdown
synthetic_patients %>%
  group_by(intubation_type) %>%
  summarise(
    n              = n(),
    fe_rate        = mean(failed_extubations > 0),
    survival_12mo  = mean(survival_12mo)
  ) %>%
  print()

# ICD chapter breakdown
synthetic_patients %>%
  count(primary_diagnosis, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = 20)


## ── 10. WRITE OUTPUTS ────────────────────────────────────────────────────────

write_csv(synthetic_patients,       "synthetic_last_extubations.csv")
write_csv(synthetic_caregiver_rate, "synthetic_caregiver_rate.csv")

cat("Written: synthetic_last_extubations.csv  (", nrow(synthetic_patients), "rows )\n")
cat("Written: synthetic_caregiver_rate.csv    (", nrow(synthetic_caregiver_rate), "rows )\n")
cat("\nThese files contain no MIMIC-IV patient data.\n")
cat("Parameters derived solely from published literature (see script header).\n")
