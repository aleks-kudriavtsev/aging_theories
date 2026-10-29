# Iteration29: three proteins beyond a conventional renal/routine core

2026-10-29. Base1e9166e03922d0666f012488ee4211b9ed0ca18c.
Recorded BEFORE new routine-domain comparisons. Positive28 is already known.
This is a new exploratory study on REUSED NHANESIII, not external/untouched
validation. Old models,28results and main remain unchanged.

## Domain and primary source definitions

First-exam NHANESIII1988–1994,40–75interview years, positive WTCYPEX6,known
eligible mortality/exam followup, original design. Same protein/oldclinical
source fingerprints as27/28, plus official1a lab/exam DAT/SAS/codebooks.
Resolve actual DAT hrefs from official index before download; no3a repeats.
EXAM bounded220MB (nominal31311*6235=195224085); LAB≤100MB. Verify byte
counts/SHA and exact INPUT ranges; FORMAT is output display, NOT an implied
inputdecimal division. Two independent readers; preserve every source record.

Common domain across ALL compared arms: original28completeclinicalcore AND
all3positivefinite reported proteins AND complete new routine inputs:
- BMI BMPBMI>0, excluding8888/blank; retain released edited/derived BMI,
  including documented source substitutions, no new imputation/cutoff. Export
  BMPWTFLG/BMPHTFLG counts, missing flags separately; do not call all measured.
- GHP HbA1c% positive finite,not8888, GHPMETH1..4. Both Diamat and affinity
  described as DCCTaligned; no duplicate diabetes-history indicator, no
  diagnostic HbA1c threshold or unverified contemporary NGSP certification.
- Original serum CEP positive finite,not8888. CDC-recalibrated
  Scr=.960*CEP−.184 must be>0. Derive adult2021 CKD-EPI creatinine-only eGFR
  using interview age/sex as in EGFR_RATIONALE.md. Estimate only, not measured
  GFR; original calibration uncertainty remains. No separate creatinine term.
- UACR=100*UBP(ug/mL)/URP(mg/dL), positivefinite UBP,not888888;
  URP>=10,not88888. ReleasedURP7.9 is a statistically-suspect/BLD replacement
  for<10mg/dL perCDC; exclude BEFORE outcomes/predictions on every arm.
  UBP0.4 below-LOD replacement is retained as RELEASED,flagged,not imputed
  anew or called aboveLOQ. No serum creatinine equation applied to urine.
No missing-value imputation, no dropping participants based on outcomes or
model advantage. Report eligible/common n,events,weighted coverage,each exclusion
reason/overlap and change versus28. No cross-domain performance comparisons.

## Prespecified core, contrasts and shapes

New binary8year weighted ridge logistic model; no frozenengine coefficients.
First6terms exactly28 age/age²,sex,current/former smoking,SBP; then
(BMI−27)/5,HbA1c,min(eGFRcr,60)/15,max(eGFRcr−60,0)/30,log2UACR.
Eleven predictor terms. eGFR is the only serum-creatinine representation,
UACR the only urinary-albumin representation. No rawcreat+eGFR duplicates,
no additional albumin/lipid/diagnosis/CHD-angina proxies or knot search.

PRIMARY: routine+ALL THREE reported protein log2values minus routine,
binary Brier8years. Secondary primary-pair AUC/calibration/bias/OE as28.
Same common cases,folds,weights for routine and routine+all3.
Recompute simple28core AND simple28core+all3 on THIS SAME common domain as
secondary descriptive comparators; never compare their absolute metrics with
saved28metrics on3857people to infer improvement or attenuation.
No single-protein winner search or panel pruning.

Predetermined exploratory assay shapes on same domain: routine+rankall3 and
routine+clipall3, same training-only ECDF/1-99%weighted log2rules as28, all
retained. Reported protein QC/LOQ/drift/storage limits from28 remain.

One limited exploratory metadata-quality sensitivity: subset of PRIMARY common
domain with BMIweightflag0 AND BMIheightflag0 AND GHPMETH1..3 AND UBP>0.4.
These conditions mean no documented anthropometric substitution, primary HPLC
and above released albumin replacement; NOT proof of per-person QC/LOQ.
Refit/evaluate both routine arms on exactly this subset with original outer
labels, no bootstrap here. Report its n/events and paired contrast; do not
compare absolute scores to main-domain scores or choose it as primary.

## Nested design, uncertainty and partition stress

8year death through96months. Stop if ANY eligible alive<96months; no silent
removal/negative outcome. Original whole-PSU primary outerhashbinary28:outer
assigned in eligible3899BEFORE complete restriction,5folds. Shared3innerfolds
per outertraining fromhash routine29:inner:{fold}, λgrid1/10/100, train-only
processing/scaling and minimalweightedinnerOOFBrier, ties1e−12→largerλ.
Same convex ridge/IRLS/finite/convergence rules as28. Normalizedtrainweights
sum rawtrainingn,includingzero-weight rows in survey bootstrap.

500paired rescaledPSUbootstrap,seed20261101,in eligibleprotein-weight3899,
restrict fixedcommon. Repeat ALL preprocessing,tuning,outerrefitting for both
primary routinearms. Save paired delta rows;2.5/97.5percentiles,count every
undefined statistic/fit. Defined-only bounds with failures are not nominal95%
coverage. Includes fitting/tuning uncertainty conditional on partition/shapes;
excludes assay/eGFRcalibration,priorselection/externaltransport uncertainty.
Calibration diagnostic fits never alter scored OOFpredictions.

LIMITED PARTITION SENSITIVITY:5 additional outcome-blind outerallocations,
SHA256('routine29:partition:{s}:'+cluster),s0..4,rankmodulo5. For each redo
nested tuning/refits for both primaryarms on same maincommon. Keep every
partition; report deltaBrier/AUC range/directions. No choosing a favorable
partition and no pooling these5dependent repeats into bootstrap intervals
or calling them5independentvalidations. No new stress option after results.

## Outputs and interpretation

Only aggregates,foldparameters,calibrationbins,500paired deltaCSV,source
manifest/metadata hashes; no participant records/IDs/predictions. Independent
percentile/hash/math audit,scientificCI and separate draftPR/evidence, no
merge/deploy/accounts/agreements/restricteddata/holdercontacts.
A positive result supports added internal reported-value information beyond
this specified routinecore; no proof of clinical validity or a universally
optimal/full panel. A null result is not biological protein uselessness.
Source primary definitions: RENAL_ASSAY_REVIEW.md,BMI_REVIEW.md,LAB_FIELDS.json,
EXAM_FIELDS.json,EGFR_RATIONALE.md. Existing28assaylimitations retained.
