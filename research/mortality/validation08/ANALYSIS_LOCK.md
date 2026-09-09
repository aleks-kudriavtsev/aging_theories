# Analysis lock: scientific validation 08

Date: 2026-09-10. This prospective analytic decision log is recorded before retrieval/inspection of participant outcomes in NHANES 2005-2008 in this project. It is not a registered clinical trial or independent protocol registration.

## Question and fixed models
Does the frozen routine laboratory block improve 10-year all-cause mortality prediction beyond clinical history, smoking, blood pressure, anthropometry, age and sex in subsequent, participant-disjoint survey cycles?

Source model: corrected 0.6_reader_corrected bundle, SHA256 7b287f3f4778b2127f2226f79d2e503e1e39087306a172e5f524a561964b07d0; code base main 6e4beccaec868b7137683f89ef912214da8883b6. No coefficients, transformations, feature cutoffs, imputation statistics or model selection will be learned from the new evaluation cycles.

Models evaluated: M0_age_sex, M0c_clinical_only_posthoc, M1_routine. M0c was post hoc in iteration04; its selection is now fixed before this evaluation. M1 includes all M0c clinical features plus total cholesterol, HDL, HbA1c, creatinine-derived eGFR, UACR, albumin, RDW, WBC and CRP. Archived extra proteins are not assumed available in the new cycles and are not silently imputed.

## Data partition
Model development: NHANES 1999-2002, corrected specimen subset (n=4549). Calibration: already examined 2003-2004 corrected specimen subset (n=2381); its previous validation results are historical, not fresh evidence. New evaluation: NHANES 2005-2006 and 2007-2008, baseline age 40-79, mortality-linkage eligible, known nonnegative examination follow-up, positive MEC examination weights and valid recorded binary sex. All eligible persons are retained with frozen training imputation; complete-case analysis is a prespecified sensitivity. Data are joined 1:1 on SEQN; no cross-partition duplicates are allowed. Both cycles must be retained regardless of results.

The target population is the eligible analytic survey sample; selection into linkage and differences from the earlier archived-specimen subset are reported. The model has historical US transport scope, not RU/DE calibration.

## Measurement harmonization
Use only CDC-prescribed cycle-specific calibration; document it before outcome evaluation. Read SAS XPORT raw zeros using xpt_checked. Reproduce original feature algebra. Creatinine, HDL and glucose/glycated haemoglobin conventions and questionnaire skip logic are checked against each cycle's public codebook. Source checksums and URLs are archived. Mortality reader is independently checked against CDC's fixed-width R specification/manual substring extraction. Zero-month event times, if any, use the previous 0.5-month convention and are reported.

## Endpoint and horizons
Primary: all-cause death within 10 years after the examination; secondary 1 and 5 years. Includes infections and external causes. It must not be labelled noninfectious mortality. Binary fixed-horizon evaluation is allowed only if every surviving participant has ascertainment through that horizon; otherwise the planned binary analysis stops and a versioned amendment is required. No 20-year extrapolation. No postmortem field is a predictor.

## Calibration strategy fixed before evaluation
Keep frozen predictions as the principal discrimination comparison. Additionally fit ONE multiplicative cumulative-hazard factor for each recorded sex and model using only 2003-2004 follow-up administratively capped at 10 years:

factor_s = sum_i(w_i * death_i_within10) / sum_i(w_i * original_cumulative_hazard_i(min(T_i,10))).

The sex-specific factors multiply all three baseline interval hazards, preserve within-sex ranks and monotone time probabilities, and do not change feature coefficients. Report before/after evaluation metrics, not just whichever is favourable. Never recalibrate on 2005-2008.

## Evaluation and uncertainty
Primary pair: M1 minus M0c fixed-horizon survey-weighted AUC at 10 years. Also report Brier difference, observed/expected, mean prediction, calibration-in-the-large, calibration slope, calibration bins, missingness and clipping rates. Pooled new-cycle weights are WTMEC2YR/2; stratum IDs include survey cycle. Paired rescaled bootstrap: 1000 replicates, resample m_h-1 PSU with replacement within each observed stratum, multiply by m_h/(m_h-1), RNG seed 20260910. Report percentile 2.5/97.5 bounds conditional on frozen models and fitted calibration factors. Do not describe these as accounting for training or calibration-estimation uncertainty.

Prespecified subgroups: female/male; age 40-59/60-79; individual survey cycles. Report event counts. Subgroup analyses are descriptive; no best-subgroup selection. Calibration regression is diagnostic only. Compare original and sex-recalibrated M1/M0c without model revision based on evaluation results.

## Scientific interpretation and release
Evidence is a participant-disjoint temporal transport test, not a randomized intervention or proof of biomarker causality. Report both favourable and unfavourable results in the technical report. A short conference abstract may concentrate on achieved results but must accurately name endpoint, population, design and estimates; no invented genetic/omics measurements. Clinical deployment and national cause-specific probabilities are not enabled by this analysis.

Code, aggregate metrics, source manifest and generated abstract facts may enter GitHub. Individual participant data/predictions remain local. Public website deployment is not needed to establish scientific validity. Any deviation must be recorded with reason and whether outcomes were inspected at that point; no overwritten lock.
