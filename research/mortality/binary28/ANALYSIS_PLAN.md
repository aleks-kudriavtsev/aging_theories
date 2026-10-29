# Iteration28: new NHANESIII eight-year binary model, reported protein values

2026-10-29. Base e86693194c6ed983838902aa9534e0d2104af4cd.
Recorded BEFORE new comparative model results. Qualification27 n3864/deaths614
is known. This is a NEW logistic model, not transport/recalibration of frozen
compact4/expanded13, and not their cause-specific piecewise-constant engine.
No old coefficients or predicted risks enter new model fitting.

## Population, variables and measurement interpretation

Same official public FIRST-exam NHANESIII1988–1994 sources/mortality2019 as27.
Age40–75 interview years, known eligible exam mortality, positive finite
WTCYPEX6, valid original design. Verify exact scientific input fingerprints
against saved27; old qualification domain3899 and joint-protein3864 must match.
Retain common participants with all3 positive finite reported proteins AND
complete interpretable clinical core, identically for every arm. No imputation.
No second-exam37a, no protected data, no new source agreement/account.

Simple clinical core: interview age HSAGEIR/HSAGEU2; HSSEX; smoking from HAR1
(>=100 lifetime cigarettes) and HAR3(current). Never=HAR1=2, current=1/1,
former=1/2; otherwise missing. HAZMNK1R overall original household+MEC mean
SBP(mmHg), up to6 readings, with HAZNOK1R in1..6;888/88/blank are missing.
Observed69–246 SBP frequency range is not a clinical cutoff. Require positive
finite SBP and valid number of readings, no frequency-range clipping.
Covariates: (age−60)/10, its square, male, current, former, (SBP−130)/20.
Age/sex/smoking/SBP is deliberately limited: no BMI, diabetes, renal function,
CVD or CHD/angina proxies. Added proteins may capture existing renal disease;
no incremental-value claim beyond a full routine/renal clinical model.
Exact positions/definitions: CLINICAL_FIELDS.json and CLINICAL_REVIEW.md.

SSCYSTAT(CST3),SSB2M,SSBTP(PTGDS), released field units mg/L; primary log2
reported values. New within-cohort fits and weighted training-only scaling,
no imported concentration coefficients or clinical cutoffs. Original Cys drift
and generic prose range-unit contradictions affect absolute interpretation;
no unverified standardization/x10/LOQ clipping is applied. Molecular identity,
specimen/visit and method are established; individual QC/LOQ flags are not.
Thus this is explicitly limited reported-measurement prognostic analysis,
not analytical/clinical validation. ASSAY_REVIEW.md records primary evidence.

## Model arms and primary contrast

Primary contrast: core+ALL THREE proteins minus core, binary Brier8years
(negative favors addition). Both arms use exactly the common domain/folds.
Secondary primary-comparison descriptors: weighted AUC, mean observed/predicted,
bias percentage points,O/E, diagnostic calibration intercept/slope.
All-cause death through96months=1, later death/alive>=96months=0.
ANY alive followup<96months stops model run: no silent removal or false survivor.
No threshold tuning, interactions, feature selection or alternative horizons.

Exploratory arms, all retained without choosing a winner: core+each single
protein (3), core+all3 weighted-ECDF proteins, core+all3 clipped log2 proteins.
ECDF Ftrain(x)=sum training weights for values<=x / sum training weights,
right-continuous, zero-weight records excluded. It is invariant to strictly
monotone common recalibration and does not establish assay validity.
Clip arm uses weighted train-only log2 1st/99th percentiles (weighted midpoint
CDF interpolation). Same domain; no record dropped as a suspected QC failure.

## Nested whole-PSU evaluation and prespecified regularization

Only masked(stratum,PSU) enters deterministic partitions. Sort all eligible
cluster keys by SHA256('binary28:outer:'+key), assign rank modulo5 BEFORE
common-case restriction. Each outer training sample forms3 inner whole-PSU
folds by SHA256('binary28:inner:{outerfold}:'+key), rank modulo3. Same folds
for all arms. No subject/PSU crosses either boundary. Folds are fixed through
bootstrap; this is internal nested CV, not independent external validation.

Weighted ridge logistic objective sum(w_normalized*(logaddexp(0,eta)−y*eta))
+lambda/2*sum(nonintercept_beta²). Weights normalized to sum=raw training n,
including zero-weight rows in resamples. Weighted train-only mean/SD for all
features (constant column SD=1, centered, explicitly flagged). Preprocessing
is refit inside every inner/outer train and bootstrap, never on heldout data.
Lambda grid[1,10,100]; choose minimal aggregate weighted inner-OOF Brier;
ties within1e−12 choose larger lambda. Refit selected lambda on outer train.
Newton/IRLS with finite objective/gradient/Hessian, line search,max100steps,
gradient tolerance1e−6 (final accepted<=1e−4 with tiny step). Failure is reported,
not rescued with synthetic events or an unplanned model.

## Uncertainty and diagnostics

500 paired rescaled-PSU bootstrap draws, seed20261031, drawn in full eligible
protein-weight domain3899 then restricted to fixed common domain. For BOTH
primary arms repeat ALL inner preprocessing,lambda selection and outer refits
with replicate weights, score using same replicate weights. Save aggregate
replicate metrics/paired deltas.2.5/97.5percentiles, report undefined/failure
counts; defined-only interval with failures is not nominal95% coverage.
Intervals include update-model fitting/tuning uncertainty conditional on fixed
partition,forms,source and chosen hypotheses. They exclude assay calibration,
prior programme selection, fold-allocation uncertainty and external transport.
500draws chosen for full nested-refit computation; endpoint MonteCarlo precision
is limited. No claim of complete individual risk uncertainty.

Calibration diagnostics fit weighted intercept-only offset and intercept+slope
against OOF logit; they do not recalibrate evaluated predictions. Report5 risk
quintile calibration bins descriptively. Bootstrap recomputes diagnostics for
primary arms; exploratory arms have point estimates/sensitivity only, no
superiority claim from their relative ranking. Report folds,support,coefficients,
preprocessing/lambda choices/convergence; no individual predictions/IDs.

## Export and decision

Only aggregate results,bootstrap rows,fold parameters,source manifest,hashes
of scientific code/plan and software versions. Public raw sources remain runner
inputs. Independent aggregate/provenance/math audit after execution.
Model remains research-only regardless of sign. Favorable internal result must
be tested with frozen complete pipeline on an independently authorized cohort;
negative result does not prove absence of biological/prognostic protein value.
No merge/deploy, no old models/results changed, new draftPRs only.

Primary sources and justified limitations are in ASSAY_REVIEW.md and
CLINICAL_REVIEW.md; all actual URLs/SHA in cohort27 SOURCE_MANIFEST.json.
