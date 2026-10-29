# Panel25: urate and wider MEC domain (2026-10-29)

Base: 4d04850231c98a0b7b07c7762d8dee9de9c2a2a5, research/nfl-robustness-joint-23.
Plan written before new individual data or comparative metrics were obtained.
Exploratory reuse of NHANES 2013–2014, NOT an untouched holdout or external validation.
Prior NfL results are known. No feature/partition/hyperparameter selection by results.

Primary extension: serum uric acid LBXSUA (mg/dL), log2 positive values,
linear coefficient with ridge=1, added to a calibration intercept multiplying
frozen expanded13 total hazard. Rationale: a directly measured renal/metabolic
analyte absent from expanded13, available in the same BIOPRO_H examination;
this is a test of incremental information, not a claim of established benefit.
No spline search, threshold search, imputation or additional analyte screening.

Expand from the NfL subsample to mortality-eligible examined adults aged 40–75
with positive finite WTMEC2YR and complete frozen predictors plus urate.
Do not require NfL measurement, do not use WTSSNH2Y for this broader domain.
Retain the entire eligible MEC domain for survey bootstrap before restricting
complete profiles. Report complete-case exclusions; no claim of unbiased
missing-data handling. Frozen baseline is the exact expanded13 JSON/SHA from
NfL22 (ten cause components, two time intervals). Older five-cause models are
different configurations; no Cox model is fitted.

First reproduce NfL22 n=1270/events=37 and its calibration/NfL OOF Brier/AUC
with the existing code on the retrieved inputs. Abort on discrepancy.
Primary horizon four years, binary all-cause death vs alive; fail closed if
any survivor is censored before four years. No multi-state Brier comparison.
Use the original outcome-blind five whole-PSU folds, assigned on the broader
eligible domain before complete-case restriction. Compare raw frozen,
calibration-only and calibration+urate on the identical complete domain.
Train using exposure to death/horizon, score using horizon cumulative hazard.
Primary metric delta Brier (urate minus calibration); secondary delta AUC.

1000 paired rescaled-PSU bootstrap draws, seed 20261029. Refit both updates
inside original folds using replicate weights; never move PSU copies between
training and scoring. Report failures and undefined statistics separately;
percentiles of defined draws with failures are not nominal 95% intervals.
Intervals condition on fixed baseline/partition and exclude baseline training,
feature-selection history and external transport uncertainty.
No model activation even if the exploratory point estimate improves.
Export only aggregates, exact source manifest, plan/code/model hashes and
software versions. Never export rows, IDs or individual predictions.

Source: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/BIOPRO_H.htm
CDC describes LBXSUA as uric acid in mg/dL, not the SI variable LBDSUASI.
Sources for all base components and mortality are the existing nfl22 retrieval
workflow; no new agreements or restricted cohorts are required.
