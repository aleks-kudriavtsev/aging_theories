# Iteration33B: separate fixed-horizon marker comparisons for cancer and heart death

Fixed before any new33B predictions, comparative metrics or model fits.
Source base daf4ebd9ec7deb18b6d8dc4e8045c3d38d869d55, stage33A protocol
82eacda12eced3c0105019fee09236b647821f45 and aggregate gate13c8d9b42f510ffee6e7b04476e545b02f9496ec08f53b585263950badf180b5.
This continues the existing research task. Only new research/evidence files are
in scope; the user's laboratory report v3 and presentation are out of scope.

## Outcomes and status of inference

Two separately fitted targets, fixed order: cancer002 (malignant neoplasms,
176 eight-year events) and heart001 (diseases of heart,148 events). Selection
follows the full-cohort outcome/burden support gate33A, not outcome-blind
selection. All cycles, samples and proteins have prior programme use. This is
exploratory internal model comparison, not independent validation, prospective
as-of deployment, confirmatory hypothesis testing or transport of frozen30.

The same6044 participants, original MEC exam origin, horizon96months inclusive,
recorded binary vital status and complete public cause ascertainment at that
horizon. For k, Yk=I(dead=1, t_exam<=96, UCOD_LEADING=k). Deaths at0/96months
count. Controls for weighted binary AUC are **all known noncases**: early known
competing deaths, deaths after96months, and survivors observed through96months.
This is not discrimination against survivors alone; it answers the specified
eight-year cumulative-incidence question. Cancer controls5868 and heart5896
are prior aggregate consistency checks. An early competing death is a known
noncase, neither censored nor considered alive. Any unknown early common cause,
early-censored survivor, invalid code/status/time or failed inherited gate
stops all dependent comparisons, without exclusion or invented zero labels.

These binary regressions estimate each published group's fixed-horizon CIF;
they are not cause-specific hazards, Cox or Fine–Gray models. Coefficients are
not hazard ratios. Do not sum independent probabilities into all-cause risk,
survival or a coherent state vector. Joint integration and Fine–Gray stay
deferred; the historical five-cause piecewise exponential ridge/CIF engine and
all-cause baselines are preserved. No exclusive marker-to-cause assignment,
causal interpretation, panel optimum, or statistically established difference
between marker effects across causes. No between-cause interaction/difference
contrast is performed; comparing signs, sizes or overlap of intervals cannot
establish one.

## Source, domain, measurement and design

Exactly the stage33A/31 primary common domain: original NHANES1999–2004 MEC,
screening age40–75, eligible public2019 mortality linkage, known status and
finite nonnegative exam followup, valid design, positive corrected SSCARD
weight, complete11 routine features and positive finite SSCYST/SSB2M values,
both detection flags0. Eligible6546, common6044, all-cause8y586, no unknown
common cause. All eligible survivors require>=96months. Source/time/cause
perturbation, assumed linkage status, missingness/complete-case selection and
historical assay calibration/QC/LOD/storage limitations remain.

Reuse the exact50 official inputs/25 pinned historical data files, corrected
SSCARDJuly2022, retained primary endpoint definitions and two-reader mortality
audit. No retry of the two prior-denied PDFs, alternate mirror, ZIP, account,
agreement or restricted data. Current source manifest and wrapper hashes are
exported. Stage33A gate must reproduce its scientific count/weight/fold/support
aggregates before any fitting; no changed-domain score comparison. Projection
keys, definitions, counts, statuses and fold fingerprints must match exactly.
Derived floating-point arithmetic is compared with predeclared relative1e-12,
absolute1e-9 tolerance solely for platform rounding; this is not permission to
alter source values, the exact native-to-W6 formula or any measurement gate.
Causes003/004/005 failed and are not fitted. Do not change horizon, thresholds,
folds, masks or group definitions to admit them.

Pooled six-year corrected subsample weights: (2/3)WTSSCB4Y for1999–2002,
(1/3)WTSSCB2Y for2003–04. Do not use native32 temporal weights or MEC weights.
The estimand is this received complete-case/measurable domain; weighting does
not undo linkage/item-missing/flag selection. Masked43 strata/87PSUs, including
stratum1's third unit; require full/eligible/common inventories to agree.

Eleven routine terms, definitions, units and fixed transformations exactly31:
age linear/quadratic,sex,current/former smoking,SBP,BMI,HbA1c,creatinine-only
eGFR low/high andlog2UACR. No duplicate raw renal predictor. Creatinine1999
1.013*raw+0.147,2001/2003 no additional adjustment; no NHANESIII equation.
Proteins ordered Cys,B2M, log2 reported serum mg/L; no new IFCC claim/correction,
assay harmonisation, imputation, threshold, clipping, interaction or tail fit.
Detection flag0 is not an individual QC/LOQ certificate.

## Arms, contrasts, metrics and multiplicity

For each cause independently, four fixed arms: routine, routine_b2m,
routine_cys, routine_cys_b2m. Every cause and arm uses identical common rows,
weights and folds. **Primary per-cause contrast:** routine_b2m minus routine,
primary metric survey-weighted8y Brier. **Secondary control:**
routine_cys_b2m minus routine_cys. AUC with the controls defined above is a
secondary diagnostic. Report all four arm scores and both contrasts, no
post-hoc winning protein/model selection or absolute Brier ranking across
different-prevalence causes.

Report observed/predicted risk, bias percentage points,O/E, calibration-in-large
(offset), unpenalized intercept+slope and unchanged-risk weighted quintile bins.
Use preserved binary28 metric definitions and diagnostic clipping only for
logits, never evaluated probabilities. Diagnostics do not recalibrate scores.
Report every undefined/nonconverged metric separately from predictive fit errors.

There are two outcomes, two contrasts and multiple diagnostics. All95% bootstrap
intervals are **marginal, exploratory, conditional percentile summaries**, not
simultaneous familywise coverage or multiplicity-adjusted confirmation. No
p-values, significance/winner declaration, favourable-cause selection or new
positive primary endpoint after results. Primary/secondary status is fixed;
all results, including null/worse ones, are retained. The bootstrap does not
include full-cohort endpoint selection, programme history, assay or transport
uncertainty. No additional interval level is chosen after the metrics.

## Nested training and paired full-refit uncertainty

Preserve the outcome-independent31 maps: five outer whole-masked-PSU folds,
hash salt b2m31:outer on all eligible units before common masking; three inner
folds within each outer training, b2m31:inner:{outer}. No fold search or
split-sensitivity variants. Every actual training split learns weighted mean
and populationSD (constant scale1). Ridge lambda1/10/100 minimizes pooled
weighted inner-heldout Brier; ties1e-12 choose larger lambda. Fresh outer
training fit after tuning, intercept unpenalized, likelihood weights normalized
to the training row count including retained zero-weight resample rows, exactly
the preserved31 numerical rule. No scoring labels enter preprocessing/tuning.

500 attempted paired rescaled within-stratum PSU bootstrap draws, RNG seed
**20261105**. Draw from the full eligible frame using the preserved31 rule:
m-1 sampled PSUs with replacement within each stratum of m units, multiplier
count*m/(m-1). The **same draw** is used for both causes and every arm; common
mask and original rows/PSU/fold identities remain fixed, including zero-weight
rows. Relearn preprocessing, lambda and outer fits on every draw. No fixed-
prediction interval, redraw, successful-only replacement, interval pooling or
seed/resample-count choice after results.

Arm/cause/metric-local failures remain in all500 attempts. AUC needs both
positive-weight classes; calibration may fail despite valid Brier/AUC. Export
valid/undefined counts and reasons; defined-only bounds do not claim full
nominal95% coverage. Percentiles2.5/97.5 use linear interpolation over finite
defined values, or null when none. All arms on a successful replicate score
unchanged OOF predictions on that same draw's weights. No extra subgroup,
assay-flag, horizon or marker variant is evaluated here.

## Reproducibility and aggregate-only export

Commit protocol/JSON contract/lock before predictions. New wrapper importer and
analyzer; all99 prior scientific modules and old models/results/main preserved.
Synthetic tests verify competing/unknown outcomes, common/draw/fold identity,
training-only nesting, failure retention and aggregate export. Scientific CI
runs source-qualified analysis and tests on Python3.11/3.13; independent audit
recomputes interval/contrast/accounting arithmetic and source/code provenance
from completed plaintext-log aggregates, disclosing no raw reprocessing.

Export RESULTS_SUMMARY.json,FOLD_PARAMETERS.json,CALIBRATION_BINS.json,
SOURCE_MANIFEST.json,bootstrap_aggregates.csv (500wide rows containing both
causes' metrics/statuses, shared weight/event-weight summaries and an aggregate
SHA256 digest of the aligned resample-weight bytes, without exporting the
weight vector or participant identities), and
paired_bootstrap.csv (2000rows:500x2causesx2contrasts). No participant IDs,
individual labels/probabilities, raw files or bootstrap coefficient vectors.
Point-fit aggregate coefficients/preprocessing are allowed research artifacts,
not a deployed/frozen replacement. Preserve exact marked stdout blocks and file
SHA, plan/source/code manifests, independent reports and commands in a separate
draft PR/evidence update. No merge,approval,deployment or changes to laboratory
reportv3/presentation.
