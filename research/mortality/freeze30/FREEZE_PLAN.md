# Iteration 30: freeze for a future external evaluation

Recorded before final-training output. This is artifact construction, not a new
prediction-performance comparison. Iteration29 results are known. No arm,
protein, transform, partition or coefficient is selected using new validation
outcomes. Old models, registries, UI and main remain unchanged.

## Training and final algorithm

Base: aging_theories 3f1a45b2bef3c77e0533d65e17707fb83f263276. Use exactly the
routine29 common domain: first-exam NHANES III, interview age 40–75, protein
subsample WTCYPEX6 >0, linked eligible mortality, all 11 routine terms and all
three positive finite reported proteins. Expected 3730 participants /563 deaths
within 96 months from examination. Preserve all routine29 source/metadata gates.
No new imputation, tail clipping, winsorization, recalibration or domain search.

Freeze exactly two arms: routine and routine_all3. Eleven routine terms, in order:
(age−60)/10, its square, male, current smoking, former smoking, (SBP−130)/20,
(BMI−27)/5, HbA1c%, min(eGFRcr,60)/15, max(eGFRcr−60,0)/30, log2(UACR).
The second arm appends log2(Cys), log2(B2M), log2(BTP). eGFRcr uses the already
reviewed 2021 race-free CKD-EPI creatinine formula on 0.960*original CEP−0.184
mg/dL, serum only. UACR=100*UBP(ug/mL)/URP(mg/dL); URP<10 remains ineligible.
Keep released UBP0.4 and official BMI edits; these are not QC/LOQ certificates.

For each arm, select ridge lambda from 1,10,100 using three whole-PSU folds on
the entire training domain, salt `freeze30:final-inner:v1`. Weighted Brier of
training-only inner heldout predictions is the tuning criterion; within 1e-12
of minimum choose the larger lambda. Each inner training split learns its own
weighted mean/population SD; constant features have SD1. Final preprocessing
and logistic ridge fit are learned afresh on all 3730 training records. Weights
are WTCYPEX6; the numerical likelihood normalizes them to the training n.
Unpenalized intercept, existing binary28 Newton/convergence rules. Same inner
folds and training participants for both arms. Export aggregate inner scores,
fold counts, coefficients, intercepts, weighted means/SDs, convergence and source
hashes. Inner tuning scores and the final-fit objective are training diagnostics,
not a new performance estimate. Export no training predictions or participant IDs.

## Scoring contract

Prediction target: probability of all-cause death by 96 months from exam, only
for the frozen research population and measurement representation. No intercept
update, refit, local scaling, imputation or tuning at scoring. Both arms require
the same complete input record, including all proteins. Inputs must name the
exact schema, contract hash, mapping and units. Unknown/additional fields, boolean
numbers, nonfinite/null/sentinel values, unknown mappings, unsupported ages and
undocumented creatinine/protein transformations are rejected before scoring.

The only initial accepted mapping is the original NHANES III first-exam reported
representation. A modern IDMS creatinine value must not be passed as historical
CEP. Modern standardized Cys, relative NPX/RFU, a second examination and other
cohorts are not silently mapped. No external mapping is approved by this freeze.
New mappings need a separately reviewed source/assay bridge, a new pinned mapping
revision and evaluation protocol committed before external outcome comparison.
Schema compliance is not authentication of provenance or legal access.

## External qualification

Review already catalogued cohorts first. Metadata-only evidence must distinguish
open person-level release, article aggregates and controlled products. Require
an independent baseline sample, exact same-visit variables, all three assays,
8-year all-cause endpoint/time origin, appropriate weights/design, access allowed
in the actual actor/environment and reviewed measurement mapping. A second visit
of training participants is not independent. Short/unknown followup cannot be
treated as complete binary status. Unknowns block; declarations cannot grant
access or validate assays. No accounts, agreements, restricted retrieval or
holder contacts. If no candidate passes, publish a concrete blocker matrix and
gate results; do not fabricate validation metrics.

## Reproducibility, tests and limits

Pin this plan, CONTRACT.json, iteration29 result/source manifest and source code.
Run final training twice in the same pinned environment and require identical
aggregate artifact bytes. Verify adapter predictions against existing numerical
primitives on temporary training data, exporting only a maximum difference and
counts. Tests cover formula/feature order, units, schema, sentinels, mapping,
tampering, whole-PSU inner separation, train-only transforms, tuning tie rule,
replay and no mutation of previous results. Independent audit will recompute
artifact hashes and predictions on synthetic fixtures, without raw records.

All prior reported-value assay/QC/LOQ, storage, historical calibration and
complete-case selection limits remain. Final artifacts are research-only;
clinical readiness and independent external validation remain false. Frozen
construction is not a claim that all proteins are useful in every cohort.
