WITH tube_events AS (

  -- Combine explicit and inferred events with a 'source' label
  SELECT *
  FROM (
    -- Explicit events
    SELECT
      p.subject_id,
      p.hadm_id,
      p.stay_id,
      p.caregiver_id,
      CAST(p.itemid AS STRING) AS itemid,
      'tube_event' AS event_type,
      p.starttime AS event_time,
      CASE
        WHEN p.itemid = 224385 THEN 'intubation'
        WHEN p.itemid = 227194 THEN 'extubation'
        WHEN p.itemid = 225468 THEN 'unplanned_extubation_patient'
        WHEN p.itemid = 225477 THEN 'unplanned_extubation_nonpatient'
      END AS tube_event,
      'explicit' AS source
    FROM `physionet-data.mimiciv_3_1_icu.procedureevents` p
    WHERE p.itemid IN (224385, 227194, 225468, 225477)

    UNION ALL

    -- Inferred intubation events
    SELECT
      i.subject_id,
      i.hadm_id,
      i.stay_id,
      NULL AS caregiver_id,
      NULL AS itemid,
      'tube_event' AS event_type,
      v.starttime AS event_time,
      'intubation' AS tube_event,
      'inferred' AS source
    FROM `mythical-legend-456217-e9.Failed_Extubation_Rate.ventilation` v
    JOIN `physionet-data.mimiciv_3_1_icu.icustays` i ON v.stay_id = i.stay_id
    WHERE v.ventilation_status = 'InvasiveVent'
      AND DATETIME_DIFF(v.endtime, v.starttime, HOUR) >= 4 --adjust to accept shorter inferred vents

    UNION ALL

    -- Inferred extubation events
    SELECT
      i.subject_id,
      i.hadm_id,
      i.stay_id,
      NULL AS caregiver_id,
      NULL AS itemid,
      'tube_event' AS event_type,
      v.endtime AS event_time,
      'extubation' AS tube_event,
      'inferred' AS source
    FROM `mythical-legend-456217-e9.Failed_Extubation_Rate.ventilation` v
    JOIN `physionet-data.mimiciv_3_1_icu.icustays` i ON v.stay_id = i.stay_id
    WHERE v.ventilation_status = 'InvasiveVent'
      AND DATETIME_DIFF(v.endtime, v.starttime, HOUR) >= 4 --adjust to accept shorter inferred vents
  )

  -- Deduplicate inferred events near explicit events of the same type
  QUALIFY NOT (
    source = 'inferred' AND (
      (LAG(tube_event) OVER w = tube_event AND LAG(source) OVER w = 'explicit' AND TIMESTAMP_DIFF(event_time, LAG(event_time) OVER w, MINUTE) BETWEEN 0 AND 30)
      OR
      (LEAD(tube_event) OVER w = tube_event AND LEAD(source) OVER w = 'explicit' AND TIMESTAMP_DIFF(LEAD(event_time) OVER w, event_time, MINUTE) BETWEEN 0 AND 30)
    )
  )

  WINDOW w AS (PARTITION BY subject_id ORDER BY event_time)

),

ranked_events AS (
  SELECT
    *,
    LAG(tube_event) OVER (PARTITION BY subject_id ORDER BY event_time) AS prev_tube_event,
    LAG(event_time) OVER (PARTITION BY subject_id ORDER BY event_time) AS prev_event_time,
    LAG(source) OVER (PARTITION BY subject_id ORDER BY event_time) AS prev_event_source,
    LEAD(event_time) OVER (PARTITION BY subject_id ORDER BY event_time) AS next_event_time,
    LEAD(tube_event) OVER (PARTITION BY subject_id ORDER BY event_time) AS next_tube_event,
    LEAD(source) OVER (PARTITION BY subject_id ORDER BY event_time) AS next_event_source
  FROM tube_events
),



extubations AS (
  SELECT
    re.subject_id,
    re.hadm_id,
    re.stay_id,
    COALESCE(re.caregiver_id, ce_near.caregiver_id, cg_fallback.caregiver_id) AS caregiver_id,
    re.event_time,
    re.tube_event,
    re.prev_event_time,
    re.prev_tube_event,
    re.next_event_time,
    re.next_tube_event,
    re.source,
    IF(re.prev_tube_event = 'intubation', 'true', 'false') AS recorded_intubation,
    IF(
      re.next_tube_event = 'intubation'
      AND TIMESTAMP_DIFF(re.next_event_time, re.event_time, HOUR) <= 72
      --AND re.source = 'explicit'
      --AND re.next_event_source = 'explicit'
      ,
      1, 0
    ) AS failed_extubation_flag,

    -- Calculate vent_hours: from prev intubation OR matched ventilation segment
    CASE
      WHEN re.prev_tube_event = 'intubation'
           AND TIMESTAMP_DIFF(re.event_time, re.prev_event_time, SECOND) > 0
        THEN TIMESTAMP_DIFF(re.event_time, re.prev_event_time, SECOND) / 3600.0
      WHEN v.starttime IS NOT NULL AND v.endtime = re.event_time
        THEN TIMESTAMP_DIFF(v.endtime, v.starttime, SECOND) / 3600.0
      ELSE NULL
    END AS vent_hours

  FROM ranked_events re
  LEFT JOIN `mythical-legend-456217-e9.Failed_Extubation_Rate.ventilation` v
    ON re.stay_id = v.stay_id
    AND v.ventilation_status = 'InvasiveVent'
    AND v.endtime = re.event_time
  
  -- If caregiver not recorded find nearest caregiver from chartevents within 15 minutes
  LEFT JOIN (SELECT
      stay_id,
      charttime,
      caregiver_id,
      ROW_NUMBER() OVER (PARTITION BY stay_id, charttime 
        ORDER BY ABS(TIMESTAMP_DIFF(charttime, charttime, MINUTE))) AS rn
    FROM `physionet-data.mimiciv_3_1_icu.chartevents`
    WHERE caregiver_id IS NOT NULL) ce_near 
    ON re.stay_id = ce_near.stay_id
    AND ABS(TIMESTAMP_DIFF(ce_near.charttime, re.event_time, MINUTE)) <= 15
    AND ce_near.rn = 1

  -- Fallback: use most frequent caregiver for that stay
  LEFT JOIN (SELECT
      stay_id,
      caregiver_id
    FROM (SELECT
        stay_id,
        caregiver_id,
        COUNT(*) AS n
      FROM `physionet-data.mimiciv_3_1_icu.chartevents`
      WHERE caregiver_id IS NOT NULL
      GROUP BY stay_id, caregiver_id
    )
  QUALIFY ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY n DESC) = 1
  ) cg_fallback 
  ON re.stay_id = cg_fallback.stay_id
  WHERE re.tube_event = 'extubation'
),

reintubation_pairs AS (
  SELECT
    e.subject_id,
    e.event_time AS extubation_time,
    MIN(i.event_time) AS reintubation_time
  FROM extubations e
  JOIN ranked_events i
    ON e.subject_id = i.subject_id
    AND i.tube_event = 'intubation'
    AND i.source = 'explicit'
    AND e.source = 'explicit'
    AND i.event_time > e.event_time
    AND TIMESTAMP_DIFF(i.event_time, e.event_time, HOUR) BETWEEN 1 AND 72
  GROUP BY e.subject_id, e.event_time
),

reintubation_summary AS (
  SELECT
    subject_id,
    COUNT(*) AS failed_extubation_count
  FROM reintubation_pairs
  GROUP BY subject_id
),


total_vent_time AS (
  SELECT
    subject_id,
    SUM(vent_hours) AS vent_hours
  FROM extubations
  GROUP BY subject_id
),

surgical_procedures AS (
  SELECT
    hadm_id,
    MIN(chartdate) AS surgery_time
  FROM `physionet-data.mimiciv_3_1_hosp.procedures_icd`
  WHERE (SUBSTR(icd_code, 1, 2) BETWEEN '01' AND '05'
       OR SUBSTR(icd_code, 1, 2) BETWEEN '06' AND '07'
       OR SUBSTR(icd_code, 1, 2) BETWEEN '36' AND '39'
       OR SUBSTR(icd_code, 1, 2) BETWEEN '42' AND '54'
       OR SUBSTR(icd_code, 1, 2) BETWEEN '76' AND '84')
  GROUP BY hadm_id
),

respiratory_diagnoses AS (
  SELECT
    hadm_id,
    COUNT(*) AS respiratory_diag_count
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (
    -- ICD-9-CM respiratory codes: 460–519
    (icd_version = 9 AND REGEXP_CONTAINS(icd_code, r'^(46[0-9]|47[0-9]|48[0-9]|49[0-9]|50[0-9]|51[0-9])'))
    
    -- ICD-10-CM respiratory codes: J00–J99
    OR (icd_version = 10 AND STARTS_WITH(icd_code, 'J'))
  )
  GROUP BY hadm_id
),

classified_intubations AS (
  SELECT
    i.subject_id,
    i.hadm_id,
    i.stay_id,
    i.event_time,
    s.surgery_time,
    r.respiratory_diag_count,
    CASE
      WHEN s.surgery_time IS NOT NULL
        AND TIMESTAMP_DIFF(TIMESTAMP(i.event_time), TIMESTAMP(s.surgery_time), HOUR)
           + COALESCE(i.vent_hours, 0) BETWEEN -12 AND 12
        THEN 'surgical'

      WHEN COALESCE(r.respiratory_diag_count, 0) > 0
        THEN 'medical-respiratory'

    ELSE 'medical-non-respiratory'
    END AS intubation_type
FROM extubations i
  LEFT JOIN surgical_procedures s ON i.hadm_id = s.hadm_id
  LEFT JOIN respiratory_diagnoses r ON i.hadm_id = r.hadm_id
),

ranked_intubations AS (
  SELECT
    i.*,
    ci.intubation_type,
    ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY i.event_time DESC) AS row_num
  FROM extubations i
  LEFT JOIN classified_intubations ci
    ON i.subject_id = ci.subject_id AND i.hadm_id = ci.hadm_id AND i.event_time = ci.event_time
  WHERE vent_hours > 0 
),

last_extubations AS (
  SELECT * EXCEPT(row_num)
  FROM ranked_intubations
  WHERE row_num = 1
),


charlson AS (
  SELECT
    subject_id,
    MAX(charlson_comorbidity_index) AS charlson_comorbidity_index
  FROM `physionet-data.mimiciv_3_1_derived.charlson`
  GROUP BY subject_id
),

norepinephrine_equivalent AS (
  SELECT
    i.subject_id,
    MAX(n.norepinephrine_equivalent_dose) AS max_norepinephrine_equivalent_dose
  FROM `physionet-data.mimiciv_3_1_derived.norepinephrine_equivalent_dose` n
  JOIN extubations i ON n.stay_id = i.stay_id
  GROUP BY i.subject_id
),

patient_info AS (
  SELECT
    subject_id,
    anchor_age,
    gender,
    dod
  FROM `physionet-data.mimiciv_3_1_hosp.patients`
),

mortality_info AS (
  SELECT
    subject_id,
    MAX(hospital_expire_flag) AS hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_hosp.admissions`
  GROUP BY subject_id
),

first_day_sofa AS (
  SELECT
    subject_id,
    MAX(sofa) AS sofa
  FROM `physionet-data.mimiciv_3_1_derived.first_day_sofa`
  GROUP BY subject_id
),

caregiver_stats AS (
  SELECT
    caregiver_id,
    COUNT(*) AS n_extubations,
    AVG(CAST(failed_extubation_flag AS FLOAT64)) AS caregiver_fe_rate,
    CASE WHEN COUNT(*) >= 30 THEN TRUE ELSE FALSE END AS is_high_volume
  FROM extubations
  WHERE caregiver_id IS NOT NULL
  GROUP BY caregiver_id
)
,
overall_fe_rate_cte AS (
  SELECT
    AVG(CAST(failed_extubation_flag AS FLOAT64)) AS overall_fe_rate
  FROM extubations
),

ccs_mappings AS (
  SELECT
  d.subject_id,
  d.hadm_id,
  d.icd_code,
  c.ccsr_category_1,
  c.ccsr_label_1
FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d
LEFT JOIN `mythical-legend-456217-e9.Failed_Extubation_Rate.ccsr_mappings` c
  ON d.icd_version = 10
 AND REPLACE(d.icd_code, '.', '') = c.icd10cm_code
),

primary_diagnoses_with_titles AS (
  SELECT
    d.hadm_id,

    -- Human-readable diagnosis (ICD-9 and ICD-10)
    MAX(
      CASE 
        WHEN dicd.icd_code IS NOT NULL AND dicd.long_title IS NOT NULL THEN 
          CONCAT(dicd.icd_code, ': ', dicd.long_title) 
        ELSE NULL 
      END
    ) AS diagnoses,

    -- CCSR mapping (ICD-10 only)
    MAX(
      CASE 
        WHEN d.icd_version = 10 AND c.ccsr_category_1 IS NOT NULL THEN 
          CONCAT(c.ccsr_category_1, ': ', c.ccsr_label_1) 
        ELSE NULL 
      END
    ) AS ccsr,

    -- ICD chapter (ICD-9 or ICD-10)
    MAX(
      CASE 
        WHEN d.icd_version = 10 THEN
          CASE
            WHEN STARTS_WITH(d.icd_code, 'A') OR STARTS_WITH(d.icd_code, 'B') THEN 'Infectious'
            WHEN STARTS_WITH(d.icd_code, 'C') OR STARTS_WITH(d.icd_code, 'D0') OR STARTS_WITH(d.icd_code, 'D1') THEN 'Neoplasms'
            WHEN STARTS_WITH(d.icd_code, 'D') THEN 'Blood/Immune'
            WHEN STARTS_WITH(d.icd_code, 'E') THEN 'Endocrine'
            WHEN STARTS_WITH(d.icd_code, 'F') THEN 'Mental'
            WHEN STARTS_WITH(d.icd_code, 'G') THEN 'Neurologic'
            WHEN STARTS_WITH(d.icd_code, 'H0') OR STARTS_WITH(d.icd_code, 'H1') THEN 'Eye'
            WHEN STARTS_WITH(d.icd_code, 'H') THEN 'Ear'
            WHEN STARTS_WITH(d.icd_code, 'I') THEN 'Circulatory'
            WHEN STARTS_WITH(d.icd_code, 'J') THEN 'Respiratory'
            WHEN STARTS_WITH(d.icd_code, 'K') THEN 'Digestive'
            WHEN STARTS_WITH(d.icd_code, 'L') THEN 'Skin'
            WHEN STARTS_WITH(d.icd_code, 'M') THEN 'Musculoskeletal'
            WHEN STARTS_WITH(d.icd_code, 'N') THEN 'Genitourinary'
            WHEN STARTS_WITH(d.icd_code, 'O') THEN 'Pregnancy'
            WHEN STARTS_WITH(d.icd_code, 'P') THEN 'Perinatal'
            WHEN STARTS_WITH(d.icd_code, 'Q') THEN 'Congenital'
            WHEN STARTS_WITH(d.icd_code, 'R') THEN 'Symptoms'
            WHEN STARTS_WITH(d.icd_code, 'S') OR STARTS_WITH(d.icd_code, 'T') THEN 'InjuryPoisoning'
            WHEN STARTS_WITH(d.icd_code, 'V') OR STARTS_WITH(d.icd_code, 'W') OR STARTS_WITH(d.icd_code, 'X') OR STARTS_WITH(d.icd_code, 'Y') THEN 'ExternalCauses'
            WHEN STARTS_WITH(d.icd_code, 'Z') THEN 'SocialFactors'
            WHEN STARTS_WITH(d.icd_code, 'U') THEN 'SpecialCodes_e.g.COVID'
            ELSE NULL
          END
        WHEN d.icd_version = 9 THEN
          CASE
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(00[1-9]|0[1-9][0-9]|1[0-1][0-9]|12[0-9]|13[0-9])') THEN 'Infectious'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(1[4-9][0-9]|2[0-3][0-9])') THEN 'Neoplasms'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(24[0-9]|25[0-9]|26[0-9]|27[0-9])') THEN 'Endocrine'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(28[0-9]|289)') THEN 'Blood/Immune'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(29[0-9]|3[0-1][0-9]|319)') THEN 'Mental'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(32[0-9]|33[0-9]|34[0-9]|35[0-9]|36[0-9]|37[0-9]|38[0-9]|389)') THEN 'Neurologic'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(39[0-9]|4[0-4][0-9]|45[0-9])') THEN 'Circulatory'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(46[0-9]|47[0-9]|48[0-9]|49[0-9]|50[0-9]|51[0-9])') THEN 'Respiratory'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(52[0-9]|53[0-9]|54[0-9]|55[0-9]|56[0-9]|57[0-9]|579)') THEN 'Digestive'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(58[0-9]|59[0-9]|60[0-9]|61[0-9]|62[0-9]|629)') THEN 'Genitourinary'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(63[0-9]|64[0-9]|65[0-9]|66[0-9]|67[0-9]|679)') THEN 'Pregnancy'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(68[0-9]|69[0-9]|70[0-9]|709)') THEN 'Skin'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(71[0-9]|72[0-9]|73[0-9]|739)') THEN 'Musculoskeletal'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(74[0-9]|75[0-9]|759)') THEN 'Congenital'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(76[0-9]|77[0-9]|779)') THEN 'Perinatal'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(78[0-9]|79[0-9]|799)') THEN 'Symptoms'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^(8[0-9][0-9]|9[0-9][0-9])') THEN 'InjuryPoisoning'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^E[0-9]{3}') THEN 'ExternalCauses'
            WHEN REGEXP_CONTAINS(d.icd_code, r'^V[0-9]{2}') THEN 'FactorsInfluencingHealth'
        ELSE 'OtherUnmapped'
    END
  ELSE NULL
      END
    ) AS icd_chapter

  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d
  LEFT JOIN `physionet-data.mimiciv_3_1_hosp.d_icd_diagnoses` dicd
    ON d.icd_code = dicd.icd_code AND d.icd_version = dicd.icd_version
  LEFT JOIN `mythical-legend-456217-e9.Failed_Extubation_Rate.ccsr_mappings` c
    ON d.icd_version = 10
   AND REPLACE(d.icd_code, '.', '') = c.icd10cm_code
  WHERE d.seq_num = 1
  GROUP BY d.hadm_id
)
,


secondary_diagnoses_with_titles AS (
  SELECT
    d.hadm_id,
    ARRAY_AGG(DISTINCT CONCAT(dicd.icd_code, ': ', dicd.long_title)) AS diagnoses
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d
  LEFT JOIN `physionet-data.mimiciv_3_1_hosp.d_icd_diagnoses` dicd
    ON d.icd_code = dicd.icd_code AND d.icd_version = dicd.icd_version
  WHERE d.seq_num = 2  
  GROUP BY d.hadm_id

),

nearest_diagnoses AS (
  SELECT
    e.subject_id,
    e.hadm_id,
    e.event_time AS extubation_time,
    e.vent_hours,
    d.icd_code,
    d.icd_version,
    dicd.long_title,
    ABS(TIMESTAMP_DIFF(e.event_time, a.admittime, HOUR)) AS time_diff_hours
  FROM last_extubations e
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON e.hadm_id = a.hadm_id
  JOIN `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d ON e.hadm_id = d.hadm_id
  LEFT JOIN `physionet-data.mimiciv_3_1_hosp.d_icd_diagnoses` dicd
    ON d.icd_code = dicd.icd_code AND d.icd_version = dicd.icd_version
),
ranked_nearest_diagnoses AS (
  SELECT *,
    ROW_NUMBER() OVER (PARTITION BY subject_id, hadm_id ORDER BY time_diff_hours) AS rn
  FROM nearest_diagnoses
),
closest_diagnosis AS (
  SELECT subject_id, hadm_id, CONCAT(icd_code, ': ', long_title) AS diagnosis_adjacent_intubation
  FROM ranked_nearest_diagnoses
  WHERE rn = 1
),



final AS (
  SELECT
    i.subject_id,
    p.anchor_age,
    p.gender,
    diag.icd_chapter AS primary_diagnosis,
    diag.diagnoses AS primary_long_diagnosis,
    diag.ccsr AS primary_ccsr,
    diag2.diagnoses AS secondary_diagnosis,
    cd.diagnosis_adjacent_intubation,
    sofa.sofa,
    c.charlson_comorbidity_index AS charlson,
    COALESCE(n.max_norepinephrine_equivalent_dose, 0) AS norepinephrine,
    va.vent_hours,
    i.event_time,
    i.intubation_type,
    i.recorded_intubation,
    i.source AS tube_event_source,
    i.caregiver_id,
    COALESCE(cs.caregiver_fe_rate, ofr.overall_fe_rate) AS caregiver_fe_rate,
    cs.n_extubations AS caregiver_n,
    COALESCE(r.failed_extubation_count, 0) AS failed_extubations,
    m.hospital_expire_flag,
    p.dod
  FROM last_extubations i
  LEFT JOIN patient_info p ON i.subject_id = p.subject_id
  LEFT JOIN reintubation_summary r 
  ON i.subject_id = r.subject_id 
  LEFT JOIN caregiver_stats cs ON i.caregiver_id = cs.caregiver_id
  LEFT JOIN charlson c ON i.subject_id = c.subject_id
  LEFT JOIN norepinephrine_equivalent n ON i.subject_id = n.subject_id
  LEFT JOIN mortality_info m ON i.subject_id = m.subject_id
  LEFT JOIN total_vent_time va ON i.subject_id = va.subject_id
  LEFT JOIN first_day_sofa sofa ON i.subject_id = sofa.subject_id
  LEFT JOIN primary_diagnoses_with_titles diag ON i.hadm_id = diag.hadm_id
  LEFT JOIN secondary_diagnoses_with_titles diag2 ON i.hadm_id = diag2.hadm_id
  LEFT JOIN closest_diagnosis cd ON i.subject_id = cd.subject_id AND i.hadm_id = cd.hadm_id
  CROSS JOIN overall_fe_rate_cte ofr 
)

SELECT *
FROM final
WHERE vent_hours IS NOT NULL 
AND vent_hours > 1 
AND vent_hours < 1000 
AND final.norepinephrine < 1
AND failed_extubations < 10