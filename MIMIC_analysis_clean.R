################################################################################
##
##  ICU Last Extubations — Analysis Script
##
##  Data source: MIMIC-IV v3.1 (BigQuery), extracted via ICU_Last_Extubations.sql
##
##  Outcome variables:
##    survival_12mo       — 12-month survival (derived from dod; TRUE if no death recorded)
##    failed_extubations  — count of reintubations within 72h of extubation
##
##  Unit of analysis:
##    Sections 1–3, 5–6: one row per patient (last extubation event)
##    Section 4:          one row per caregiver
##    Section 7:          one row per caregiver (wide diagnosis case mix)
##    Section 8:          one row per extubation event with RSBI recorded
##
##  SQL pipeline summary (ICU_Last_Extubations.sql):
##    - Combines explicit procedure events (itemids 224385, 227194, 225468, 225477)
##      with inferred events from ventilation segments (>= 4h InvasiveVent)
##    - Deduplicates inferred events within 30 min of explicit events
##    - Classifies intubation type: surgical / medical-respiratory / medical-non-respiratory
##    - Selects last extubation per patient
##    - Joins: Charlson comorbidity index, norepinephrine equivalent dose,
##             SOFA score (first day), hospital mortality, ICD chapter, CCSR category
##    - Assigns caregiver via: recorded ID → nearest chartevents (±15 min) → mode fallback
##    - Computes caregiver_fe_rate, caregiver_n at query time
##    - Final filters: vent_hours in (1, 1000), norepinephrine < 1, failed_extubations < 10
##
################################################################################


## ── 0. LIBRARIES AND SETUP ───────────────────────────────────────────────────

library(tidyverse)
library(lubridate)
library(splines)
library(glmnet)
library(randomForest)
library(caret)

set.seed(237)

# setwd("~/Projects/Extubate/r")   # adjust as needed


## ── 1. DATA INGESTION ────────────────────────────────────────────────────────
##
##  All CSVs are outputs of ICU_Last_Extubations.sql run on BigQuery.
##  tube_events and patient_timelines are intermediate outputs used only
##  for timeline visualization (Section 9).

last_extubations <- read_csv(
  "../data/last_extubations.csv",
  col_types = cols(
    subject_id          = col_character(),
    gender              = col_factor(),
    intubation_type     = col_factor(),
    caregiver_id        = col_character(),
    event_time          = col_datetime(),
    tube_event_source   = col_factor(),
    dod                 = col_datetime()
  )
)

last_extubations$survival_12mo <- is.na(last_extubations$dod)

summary(last_extubations)


## ── 2. COHORT CONSTRUCTION ───────────────────────────────────────────────────
##
##  cleandata: primary analytical cohort
##    - caregiver_n > 10: ensures caregiver_fe_rate is estimated from meaningful volume
##    - drops rows with missing primary_diagnosis or intubation_type
##    - centers continuous predictors (for interpretability of main effects in
##      models with polynomial caregiver_fe_rate terms)
##
##  supercleandata: restricted to explicit tube events only
##    - used for sensitivity checks; smaller N but higher data quality

cleandata <- last_extubations %>%
  filter(
    caregiver_n > 10,
    !is.na(primary_diagnosis),
    !is.na(intubation_type)
  ) %>%
  mutate(across(
    c(vent_hours, charlson, norepinephrine, sofa),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_centered"
  ))

cat("Cohort N:", length(unique(cleandata$subject_id)), "\n")
cat("Failed extubation rate:", mean(cleandata$failed_extubations > 0), "\n")
cat("12-month survival rate:", mean(cleandata$survival_12mo, na.rm = TRUE), "\n")

# Intubation type breakdown
cleandata %>%
  group_by(intubation_type) %>%
  summarise(
    n = n(),
    fe_rate = mean(failed_extubations > 0),
    survival_12mo = mean(survival_12mo, na.rm = TRUE)
  )

# Sensitivity cohort: explicit events only
supercleandata <- cleandata %>%
  filter(tube_event_source == "explicit")

cat("Explicit-only cohort N:", nrow(supercleandata), "\n")
cat("Explicit-only FE rate:", mean(supercleandata$failed_extubations > 0), "\n")

supercleandata %>%
  group_by(intubation_type) %>%
  summarise(fe_rate = mean(failed_extubations > 0))


## ── 3. PATIENT-LEVEL LOGISTIC REGRESSION: 12-MONTH SURVIVAL ─────────────────
##
##  Outcome: survival_12mo (TRUE = alive at 12 months)
##  Core predictors: intubation_type, vent_hours, charlson, norepinephrine, sofa
##  Key covariate: caregiver_fe_rate (caregiver's historical failed extubation rate)
##
##  Model variants:
##    logreg1    — centered predictors, poly(caregiver_fe_rate, 2)
##    logreg2    — suppressed intercept, natural spline on caregiver_fe_rate
##    logreg3    — full interaction: all predictors × poly(caregiver_fe_rate, 2)
##
##  Note: logreg1_nocenter is included to verify centering does not change
##        fitted values (coefficients change but predictions are identical).

logreg1 <- glm(
  survival_12mo ~ intubation_type +
    vent_hours_centered +
    charlson_centered +
    norepinephrine_centered +
    sofa_centered +
    poly(caregiver_fe_rate, 2),
  data = cleandata, family = binomial
)
summary(logreg1)
round(exp(logreg1$coefficients), 3)
cleandata$logreg1_fitted <- logreg1$fitted.values

# Verify: centered vs. uncentered models give same fitted values
logreg1_nocenter <- glm(
  survival_12mo ~ intubation_type +
    vent_hours +
    charlson +
    norepinephrine +
    sofa +
    poly(caregiver_fe_rate, 2),
  data = cleandata, family = binomial
)
cleandata$logreg1_nocenter_fitted <- logreg1_nocenter$fitted.values
stopifnot(max(abs(cleandata$logreg1_fitted - cleandata$logreg1_nocenter_fitted)) < 1e-8)

# Alternative: no-intercept model with natural spline on caregiver_fe_rate
logreg2 <- glm(
  survival_12mo ~ 0 + intubation_type +
    vent_hours_centered +
    charlson_centered +
    norepinephrine_centered +
    sofa_centered +
    ns(caregiver_fe_rate, 3),
  data = cleandata, family = binomial
)
summary(logreg2)
cleandata$logreg2_fitted <- logreg2$fitted.values

# Compare logreg1 vs logreg2 fitted values
ggplot(cleandata, aes(logreg1_fitted, logreg2_fitted)) +
  geom_point(alpha = 0.3) +
  geom_abline(linetype = "dashed") +
  labs(title = "Fitted value comparison: poly vs. natural spline on caregiver_fe_rate")

# Full interaction model
logreg3 <- glm(
  survival_12mo ~ 0 + (intubation_type +
    vent_hours +
    charlson +
    norepinephrine +
    sofa) * poly(caregiver_fe_rate, 2),
  data = cleandata, family = binomial
)
summary(logreg3)
cleandata$logreg3_fitted <- logreg3$fitted.values
round(exp(logreg3$coefficients), 3)


## ── 4. CAREGIVER-LEVEL ANALYSIS ──────────────────────────────────────────────
##
##  Aggregates to caregiver level to examine:
##    - Does caregiver volume predict FE rate? (volume-outcome relationship)
##    - Does caregiver FE rate predict patient survival / hospital mortality?
##
##  caregiver_diagnoses: caregiver × diagnosis breakdown (for LASSO in Section 5)
##  caregiver_rate:      one row per caregiver with summary rates

caregiver_rate <- cleandata %>%
  group_by(caregiver_id) %>%
  summarise(
    caregiver_fe_rate      = max(caregiver_fe_rate),
    caregiver_n            = max(caregiver_n),
    hospital_expire_rate   = mean(hospital_expire_flag),
    survival_12mo_rate     = mean(survival_12mo)
  )
summary(caregiver_rate)

mean_cg_fe  <- mean(caregiver_rate$caregiver_fe_rate)
med_cg_fe   <- median(caregiver_rate$caregiver_fe_rate)

# Volume vs. FE rate
ggplot(caregiver_rate, aes(caregiver_n, caregiver_fe_rate)) +
  geom_point() +
  geom_hline(yintercept = mean_cg_fe) +
  geom_hline(yintercept = med_cg_fe, linetype = "dashed") +
  geom_smooth() +
  labs(
    title = "Caregiver Volume vs. Failed Extubation Rate",
    x = "Number of extubations (caregiver_n)",
    y = "Failed extubation rate"
  )

# Volume vs. 12-month survival rate
mean_cg_surv <- mean(caregiver_rate$survival_12mo_rate)
med_cg_surv  <- median(caregiver_rate$survival_12mo_rate)

ggplot(caregiver_rate, aes(caregiver_n, survival_12mo_rate)) +
  geom_point() +
  geom_hline(yintercept = mean_cg_surv) +
  geom_hline(yintercept = med_cg_surv, linetype = "dashed") +
  geom_smooth() +
  labs(
    title = "Caregiver Volume vs. 12-Month Patient Survival Rate",
    x = "Number of extubations (caregiver_n)",
    y = "12-month survival rate"
  )

# Volume vs. hospital mortality
mean_hosp_exp <- mean(caregiver_rate$hospital_expire_rate)
med_hosp_exp  <- median(caregiver_rate$hospital_expire_rate)

ggplot(caregiver_rate, aes(caregiver_n, hospital_expire_rate)) +
  geom_point() +
  geom_hline(yintercept = mean_hosp_exp) +
  geom_hline(yintercept = med_hosp_exp, linetype = "dashed") +
  labs(
    title = "Caregiver Volume vs. Hospital Mortality Rate",
    x = "Number of extubations (caregiver_n)",
    y = "Hospital expire rate"
  )

# Caregiver × primary diagnosis breakdown (input to Section 5 LASSO)
caregiver_diagnoses <- cleandata %>%
  group_by(caregiver_id, primary_diagnosis) %>%
  summarise(
    caregiver_fe_rate      = max(caregiver_fe_rate),
    caregiver_n            = max(caregiver_n),
    diag_rate              = n() / max(caregiver_n),
    hospital_expire_rate   = mean(hospital_expire_flag),
    survival_12mo_rate     = mean(survival_12mo),
    .groups = "drop"
  ) %>%
  filter(!(primary_diagnosis %in% c("Ear")))  # too few cases to estimate reliably


## ── 5. DIAGNOSIS CASE MIX: LASSO FEATURE SELECTION ──────────────────────────
##
##  Question: which ICD chapters in a caregiver's case mix predict their FE rate?
##
##  Approach:
##    1. Pivot caregiver_diagnoses to wide format (one column per ICD chapter)
##    2. Remove near-zero-variance diagnosis columns (caret::nearZeroVar)
##    3. Fit LASSO (cv.glmnet) to select chapters predictive of caregiver_fe_rate
##    4. Manually add selected ICD chapter indicators back to cleandata
##       (positive predictors = higher FE risk; negative = lower FE risk)

wide_data <- caregiver_diagnoses %>%
  pivot_wider(
    id_cols    = c(caregiver_id, caregiver_fe_rate),
    names_from  = primary_diagnosis,
    values_from = diag_rate,
    values_fill = list(diag_rate = 0)
  )

y <- wide_data$caregiver_fe_rate
X <- wide_data %>%
  select(-caregiver_id, -caregiver_fe_rate) %>%
  as.matrix()

set.seed(123)
cvfit      <- cv.glmnet(X, y, alpha = 1)
plot(cvfit)

best_lambda <- cvfit$lambda.min
lasso_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)

coef_df <- as.data.frame(as.matrix(coef(lasso_model))) %>%
  tibble::rownames_to_column("diag_code") %>%
  filter(s0 != 0)
print(coef_df)

## ICD chapter indicators derived from LASSO output
##   Positive predictors of FE rate (higher caregiver FE risk)
cleandata$icd_digestive   <- cleandata$primary_diagnosis == "Digestive"
cleandata$icd_special     <- cleandata$primary_diagnosis == "SpecialCodes_e.g.COVID"
cleandata$icd_infectious  <- cleandata$primary_diagnosis == "Infectious"
cleandata$icd_neurologic  <- cleandata$primary_diagnosis == "Neurologic"

##   Negative predictors of FE rate (lower caregiver FE risk)
cleandata$icd_circulatory    <- cleandata$primary_diagnosis == "Circulatory"
cleandata$icd_respiratory    <- cleandata$primary_diagnosis == "Respiratory"
cleandata$icd_genitourinary  <- cleandata$primary_diagnosis == "Genitourinary"
cleandata$icd_pregnancy      <- cleandata$primary_diagnosis == "Pregnancy"
cleandata$icd_symptoms       <- cleandata$primary_diagnosis == "Symptoms"
cleandata$icd_musculoskeletal <- cleandata$primary_diagnosis == "Musculoskeletal"

## ICD chapter grouping for visualization
top_n       <- 8
top_chapter <- cleandata %>%
  count(primary_diagnosis, sort = TRUE) %>%
  slice_head(n = top_n) %>%
  pull(primary_diagnosis)

cleandata <- cleandata %>%
  mutate(
    primary_icd_grouped = case_when(
      primary_diagnosis == ""                   ~ "Missing",
      primary_diagnosis %in% top_chapter        ~ primary_diagnosis,
      TRUE                                      ~ "Other"
    ),
    primary_icd_grouped = fct_infreq(primary_icd_grouped)
  )


## ── 6. FULL MODEL WITH ICD INTERACTIONS + 5-FOLD CROSS-VALIDATION ────────────
##
##  logreg_icd: survival_12mo modeled with ICD chapter indicators interacting
##  with poly(caregiver_fe_rate, 2). This tests whether the effect of caregiver
##  FE rate on patient survival varies by diagnosis category.
##
##  5-fold CV reports: accuracy, PPV, NPV, sensitivity, specificity
##  Coefficient stability across folds assessed via SD/mean ratio.

icd_formula <- survival_12mo ~
  poly(caregiver_fe_rate, 2) * (
    intubation_type +
    vent_hours +
    charlson_centered +
    norepinephrine_centered +
    sofa_centered +
    icd_digestive +
    icd_genitourinary +
    icd_respiratory +
    icd_neurologic +
    icd_special +
    icd_infectious +
    icd_circulatory +
    icd_pregnancy +
    icd_musculoskeletal
  )

logreg_icd <- glm(icd_formula, data = cleandata, family = binomial)
summary(logreg_icd)
cleandata$logreg_icd_fitted <- predict(logreg_icd, type = "response")

# Diagnostic plots: fitted values by key predictors, faceted by ICD chapter
plot_vars <- c("caregiver_fe_rate", "charlson", "sofa", "vent_hours", "anchor_age")

for (v in plot_vars) {
  p <- ggplot(cleandata, aes_string(v, "logreg_icd_fitted",
                                     color = "survival_12mo",
                                     shape = "intubation_type")) +
    geom_point(alpha = 0.3) +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~primary_icd_grouped) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1)),
      shape = guide_legend(override.aes = list(alpha = 1))
    ) +
    labs(title = paste("Fitted survival probability vs.", v))
  print(p)
}

# Caregiver FE rate binned boxplot
cleandata$caregiver_fe_binned <- cut(
  cleandata$caregiver_fe_rate,
  breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 1),
  right  = FALSE
)
cleandata <- cleandata %>% filter(!is.na(caregiver_fe_binned))

ggplot(cleandata, aes(caregiver_fe_binned, logreg_icd_fitted, fill = survival_12mo)) +
  geom_boxplot() +
  labs(
    title = "Model-fitted survival by caregiver FE rate bin",
    x     = "Caregiver failed extubation rate (binned)",
    y     = "Predicted 12-month survival probability"
  )

# Low vs. high FE caregiver: top diagnosis comparison
low_fe  <- cleandata %>% filter(caregiver_fe_rate < 0.01)
high_fe <- cleandata %>% filter(caregiver_fe_rate > 0.09)

cat("Low FE caregiver — top CCSR diagnoses:\n")
head(sort(table(low_fe$primary_ccsr), decreasing = TRUE))

cat("High FE caregiver — top CCSR diagnoses:\n")
head(sort(table(high_fe$primary_ccsr), decreasing = TRUE))

# 5-fold cross-validation
set.seed(123)
k     <- 5
n     <- nrow(cleandata)
folds <- sample(rep(1:k, length.out = n))

accuracies   <- numeric(k)
ppvs         <- numeric(k)
npvs         <- numeric(k)
sensitivities <- numeric(k)
specificities <- numeric(k)
coefs_list   <- vector("list", k)

for (i in 1:k) {
  train_data <- cleandata[folds != i, ]
  test_data  <- cleandata[folds == i, ]

  model <- glm(icd_formula, data = train_data, family = binomial)
  coefs_list[[i]] <- coef(model)

  pred_probs <- predict(model, newdata = test_data, type = "response")
  preds      <- pred_probs > 0.5
  actuals    <- test_data$survival_12mo

  TP <- sum(preds  &  actuals, na.rm = TRUE)
  TN <- sum(!preds & !actuals, na.rm = TRUE)
  FP <- sum(preds  & !actuals, na.rm = TRUE)
  FN <- sum(!preds &  actuals, na.rm = TRUE)

  accuracies[i]    <- mean(preds == actuals, na.rm = TRUE)
  ppvs[i]          <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  npvs[i]          <- if ((TN + FN) > 0) TN / (TN + FN) else NA
  sensitivities[i] <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  specificities[i] <- if ((TN + FP) > 0) TN / (TN + FP) else NA
}

cat("Cross-validated accuracies:      ", round(accuracies, 3), "\n")
cat("Mean accuracy:                   ", round(mean(accuracies, na.rm = TRUE), 3), "\n")
cat("SD accuracy:                     ", round(sd(accuracies, na.rm = TRUE), 3), "\n")
cat("Mean PPV:                        ", round(mean(ppvs, na.rm = TRUE), 3), "\n")
cat("SD PPV:                          ", round(sd(ppvs, na.rm = TRUE), 3), "\n")
cat("Mean NPV:                        ", round(mean(npvs, na.rm = TRUE), 3), "\n")
cat("SD NPV:                          ", round(sd(npvs, na.rm = TRUE), 3), "\n")
cat("Mean Sensitivity (Recall):       ", round(mean(sensitivities, na.rm = TRUE), 3), "\n")
cat("SD Sensitivity:                  ", round(sd(sensitivities, na.rm = TRUE), 3), "\n")
cat("Mean Specificity:                ", round(mean(specificities, na.rm = TRUE), 3), "\n")
cat("SD Specificity:                  ", round(sd(specificities, na.rm = TRUE), 3), "\n")

# Coefficient stability across folds
coef_matrix  <- do.call(rbind, coefs_list)
coef_summary <- data.frame(
  term = colnames(coef_matrix),
  mean = colMeans(coef_matrix, na.rm = TRUE),
  sd   = apply(coef_matrix, 2, sd, na.rm = TRUE)
)
cat("Most stable coefficients (lowest CV):\n")
head(coef_summary[order(abs(coef_summary$sd / coef_summary$mean)), ], 10)


## ── 7. RANDOM FOREST: CAREGIVER CASE MIX → FE RATE ──────────────────────────
##
##  Unit of analysis: caregiver
##  Features: proportion of each primary diagnosis in caregiver's caseload
##  Target:   caregiver_fe_rate
##
##  Near-zero-variance diagnosis columns removed before fitting.
##  Variable importance assessed via %IncMSE and IncNodePurity.

caregiver_wide <- caregiver_diagnoses %>%
  select(caregiver_id, primary_diagnosis, diag_rate, caregiver_fe_rate) %>%
  pivot_wider(
    names_from  = primary_diagnosis,
    values_from = diag_rate,
    values_fill = list(diag_rate = 0)
  ) %>%
  group_by(caregiver_id) %>%
  summarise(across(everything(), max))

nzv          <- nearZeroVar(caregiver_wide[, -1], saveMetrics = TRUE)
filtered_vars <- caregiver_wide[, c(TRUE, !nzv$nzv)]

data_rf <- filtered_vars %>%
  select(-caregiver_id) %>%
  rename_with(~ gsub("\\[|\\]|\\/| |,", "", .x))

set.seed(123)
rf_model <- randomForest(
  caregiver_fe_rate ~ .,
  data       = data_rf,
  importance = TRUE,
  ntree      = 50
)
summary(rf_model)

importance_df <- importance(rf_model) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("diag_code") %>%
  arrange(desc(`%IncMSE`))

top_vars  <- importance_df %>% slice_max(`%IncMSE`, n = 20)
top_vars2 <- importance_df %>% slice_max(IncNodePurity, n = 20)

ggplot(top_vars, aes(x = reorder(diag_code, `%IncMSE`), y = `%IncMSE`)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 20 Diagnosis Codes Predicting Caregiver FE Rate (% Inc MSE)",
    x = "Diagnosis Code", y = "% Increase in MSE"
  )

ggplot(top_vars2, aes(x = reorder(diag_code, IncNodePurity), y = IncNodePurity)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 20 Diagnosis Codes Predicting Caregiver FE Rate (Node Purity)",
    x = "Diagnosis Code", y = "Increase in Node Purity"
  )


## ── 8. RSBI-BASED MODEL (SEPARATE COHORT) ────────────────────────────────────
##
##  Simpler logistic regression using pre-extubation vitals including RSBI
##  (Rapid Shallow Breathing Index). This cohort is separate from cleandata
##  as it requires a chart-level vitals join (extubation_vitals) and the
##  diagnoses table structured with binary flags.
##
##  Note: extubation_df is loaded from a pre-joined CSV exported from BigQuery.
##        The diagnoses column parsing below handles the array-as-string format
##        that BigQuery exports for repeated fields.

extubation_df <- read.csv("extubation.csv")

extubation_df <- extubation_df %>%
  mutate(
    diagnoses = diagnoses %>%
      str_remove_all("\\[|\\]") %>%
      str_split(",(?=[^,]*:)") %>%
      map(str_trim)
  )

summary(extubation_df)

# Model 1: diagnosis flags only
logreg_rsbi1 <- glm(
  (failed_extubation_flag == "true") ~
    has_respiratory_failure +
    has_ards +
    has_copd +
    has_pneumonia +
    has_sepsis,
  data = extubation_df, family = binomial
)
summary(logreg_rsbi1)

# Model 2: add RSBI
logreg_rsbi2 <- glm(
  (failed_extubation_flag == "true") ~
    has_respiratory_failure +
    has_copd +
    has_pneumonia +
    has_ards +
    has_sepsis +
    rsbi_final,
  data = extubation_df, family = binomial
)
summary(logreg_rsbi2)


## ── 9. APPENDIX: PATIENT TIMELINE VISUALIZATION ──────────────────────────────
##
##  Exploratory visualization of intubation/extubation event sequences.
##  Not part of primary analysis; retained for data quality checking.

tube_events <- read_csv(
  "../data/tube_events.csv",
  col_types = cols(
    event_type  = col_factor(),
    tube_event  = col_factor(),
    source      = col_factor()
  )
)

patient_timelines <- read_csv(
  "../data/patient_timelines.csv",
  col_types = cols(
    subject_id = col_character(),
    event_time = col_datetime(),
    stay_id    = col_character()
  )
)

patient_timelines <- patient_timelines %>%
  group_by(subject_id) %>%
  mutate(
    first_intubation_time         = min(event_time[tube_event == "intubation"], na.rm = TRUE),
    hours_since_first_intubation  = as.numeric(difftime(event_time, first_intubation_time, units = "hours"))
  ) %>%
  ungroup() %>%
  filter(is.finite(first_intubation_time)) %>%
  filter(event_time >= first_intubation_time)

sample_ids       <- sample(unique(patient_timelines$subject_id), 5)
sample_timelines <- patient_timelines %>% filter(subject_id %in% sample_ids)

ggplot(sample_timelines, aes(x = hours_since_first_intubation, y = factor(subject_id))) +
  geom_point(aes(color = tube_event, shape = source), size = 3) +
  labs(
    x     = "Hours since first intubation",
    y     = "Patient (subject_id)",
    title = "Timeline of Intubation and Extubation Events",
    color = "Tube Event",
    shape = "Source"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), legend.position = "bottom")
