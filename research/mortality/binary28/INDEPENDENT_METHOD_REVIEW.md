# Independent scientific/code review28 before comparative metrics

Reviewed2026-10-29: committed plan d81972a, CLINICAL_FIELDS/REVIEW,
ASSAY_REVIEW, LOCK, cohort.py, analyze.py, metrics.py and the now-present model.py
under research/mortality/binary28. Read-only inspection; no participant data,
model execution, new empirical calculation, ZIP or account action. Files were
being completed concurrently, so findings apply to the reviewed implementation,
not an unobserved later revision.

## Verdict

No critical scientific inconsistency or leakage found in the stated design and
current implementation. It is a NEW limited-core weighted binary logistic model
of eight-year all-cause mortality, not frozen compact4/expanded13 transport or
cause-specific modelling. A reported-values internal prognostic comparison is
interpretable under the explicitly retained assay, missingness, sampling and
validation limitations. A favorable result would not establish clinical assay
validity, external validation or incremental value beyond a full renal/clinical
baseline.

## Source and clinical handling

- The new reader correctly requires support through column3341, because smoking
  and released SBP/count are far beyond qualification27's first45columns.
  It isolates SAS INPUT rather than LENGTH, checks exact ranges, respects
  documented MISSOVER only after all selected fields are present, compares two
  readers exactly and retains the official20050adult records. It does not pad
  absent selected fields or infer their contents.
- All non-index scientific source bytes must match saved27 fingerprints. Plan,
  clinical/assay metadata and source-manifest lock checks precede processing.
  Reusing the27public-source adapter does not transfer old model coefficients:
  binary28 defines its own plan and model; frozen model application remains out
  of scope.
- Interview age and sex semantics are retained. HAR1/HAR3 produce the specified
  >=100-cigarette-based never/current/former categories; unknown codes are not
  made never. This does not measure every tobacco product or exposure history.
- HAZMNK1R is explicitly the original household+MEC released SBP summary, not a
  later repeat-exam measurement.888, invalid/missing reading count and blank
  remain missing. No frequency-range-based physiological clipping is introduced.
- The limited age/sex/smoking/SBP core deliberately excludes BMI, diabetes and
  renal predictors. Protein gain could largely capture existing renal disease;
  the plan and result scope make that limitation explicit.
- Common complete cases include allthree proteins for EVERY comparison, including
  single-protein exploratory arms. No imputation or model-specific population
  change was found. Alive followup<96months causes a stop before modelling;
  later deaths are known nonevents at8years; deaths<=96months are events.

## Reported protein interpretation

Primary log2 original reported values plus training-only centering/scaling is a
within-cohort prediction question, not use of a transported concentration cutoff.
Individual QC/LOQ remains unknown. No x10 repair, guessed assay range screen or
unverified cystatin standardization is applied. Fixed ECDF and clip sensitivities
are appropriately labelled exploratory and do not replace the primary result.

ECDF is computed on training positive-weight support, with a defined right-
continuous tie rule and zero below minimum/one above maximum. On log values it
has the same ordering as original concentrations and is invariant to a common
strictly increasing recalibration. It cannot fix specimen-specific or batch-
specific errors, freeze/thaw-induced reorderings, selection or confounding.
Clipped-log bounds use the prescribed weighted midpoint1st/99th percentiles;
they are robustness quantiles, not LOD/LOQ or clinical reference limits.

## Nested design and refitting uncertainty

- Outer partitions depend only on eligible PSU identity, before common restriction.
  Inner partitions use only the outer training PSU identities. Source rows and
  PSU membership remain present with zero replicate weight, so original fold
  boundaries stay fixed through bootstrap. Both levels explicitly reject overlap.
- Training-only transformation/scaling is refit within each inner split and final
  outer train. Its cache is shared across lambda values only within that same
  training split; this is safe because preprocessing is outcome-independent.
- Lambda1/10/100 is selected by aggregated weighted inner heldout Brier, with
  prescribed ties favoring stronger shrinkage. Outer outcomes do not affect
  their own model fitting or tuning. No full-sample rank or clipping fit found.
- Ridge objective, unpenalized intercept, weight normalization to raw training n
  including zero-weight rows and finite gradient/Hessian checks match the plan.
  Numerical failures are explicit rather than silently replaced.
- Bootstrap weights are drawn in the entire eligible3899domain, THEN restricted
  to fixed common cases. Both primary arms rerun preprocessing, inner tuning and
  outer training, and are scored with identical replicate weights. The paired
  interval therefore includes fitted-model/tuning variation conditional on the
  fixed fold allocation and chosen forms; it is not a fixed-prediction bootstrap.
- The interval still excludes assay calibration/measurement uncertainty, prior
  research/hypothesis selection, partition allocation and transport.500draws give
  limited percentile endpoint precision; undefined draws and diagnostics need
  explicit counts, which the exporter provides.

## Metrics and export

Brier, AUC and mean calibration all use survey/replicate weights. Diagnostic
calibration intercept/slope are fitted after OOF prediction but never alter the
scored risks. Probability clipping applies only to diagnostic logits. Single-
class, flat/separated or unconverged calibration is labelled undefined rather
than inventing a slope. Quintile bins are descriptive, tied scores stay together.

Only aggregate metrics/parameters/support and code/source provenance are emitted.
The model's internal ECDF support arrays are not exported in fold parameters;
this avoids publishing a list of participant measurement values. Model-family,
no-legacy-coefficient, internal-validation and reported-value/QC limits are explicit.

## Nonblocking auditability improvements

1. In analyze.py bootstrap_penalty_selection_counts is incremented after each arm,
   before the complete core/all3 pair succeeds. If all3 later fails, successful
   core fits from an undefined replicate remain counted. Either accumulate both
   selections only after a successful pair, or label counts as all successful
   intermediate fits, including failed-pair replicates. Metrics/interval values
   themselves are not affected by this reporting issue.
2. Fold records identify clipped-log transformation but omit its six aggregate
   clip bounds. Saving those scalar bounds would improve sensitivity auditability
   without publishing participant arrays. Full ECDF support should remain private.
3. calibration_bins.csv is artifact-only. Under the known ZIP restriction, marked
   aggregate JSON stdout for bins would permit independent audit without another
   model calculation or download workaround. Current primary scores are available
   in marked JSON stdout, so this is not a primary-result blocker.

No implementation edit was made by this reviewer. These are pre-metric findings;
no choice of transform or model arm was informed by observed comparison results.

## Follow-up read-only check on scientific code e2171f9

Before reading any empirical comparison results, inspected the resulting code:
penalty-selection counters are now updated only after both models and paired
scoring succeed; calibration bins are printed in a separate marked aggregate
JSON block. Scalar clip bounds are exported via protein_transform_parameters, while rank
support values remain private. All three original auditability items were
communicated and resolved before reading metrics. No statistical
protocol or primary contrast was changed in response to outcomes.
