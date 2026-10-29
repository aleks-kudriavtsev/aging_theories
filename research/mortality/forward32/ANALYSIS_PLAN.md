# Iteration32: fixed forward-calendar stress test of the B2M contrasts

Record before new temporal model metrics. Earlier31 results are known: B2M adds
information beyond routine, but its gain after Cys is uncertain. Both calendar
parts have prior programme use, including old M4 and31. This is exploratory
strict-calendar evaluation on reused persons, not an untouched/independent
holdout, confirmation of a winning panel or external validation of frozen30.
Preserve all earlier models, protocols, results, main and UI.

## Single fixed split and outcome

Train ONLY baseline1999–2000 plus2001–2002; score ONLY baseline2003–2004.
No reverse training, rolling split search, selection of a favorable cycle,
late refit, intercept update or recalibration. Biological baseline is exam,
not serum assay2018–2020. Eight-year all-cause outcome from2019 public LMF:
dead=1 and exam followup<=96months. Missing status remains missing; survivors
must have>=96months. No synthetic replacement of outcomes or restricted2022.

Reuse exactly b2m31 source/measurement eligibility and primary common mask:
screening age40–75, MEC examined, linked eligible, finite routine11 predictors,
positive finite reported Cys/B2M and both flags0. Same common domain in all four
arms within each calendar part. Do not add31's LLOD sensitivity or new subgroups.
Known31 common counts4020earlier/2024later are historical planning information;
new temporal event counts and comparative metrics have not been inspected.

## Weights and target

Use original corrected SSCARD_A **WTSSCB4Y** for pooled1999–2002 training;
use original corrected **WTSSCB2Y** for2003–2004 evaluation. Explicitly replace
the old frame's W6 `weight`; retain raw fields/W6 only for provenance. Do not
drop late cases because their irrelevant WTSSCB4Y is missing. Do not reuse
W6 labels, normalize each early cycle separately or pool train/test into one
score. Survey weights in scoring remain raw native period weights; only the
training ridge likelihood is normalized to training n as already specified.

The empirical estimand is the received complete-case/measurable-protein domain
of each respective period. Native weights support the source subsample, not
an automatic correction for extra linkage/missing/LOD selection. No unqualified
national-population validation claim, missingness reweighting or transport to
2026/other countries. Report eligible/common native weight sums, raw W4/W2
support, W6 comparison, coverage and Kish effective n separately for each part.

## Source/design/support gate before model fits

Pin90 prior scientific Python hashes,31 metadata/results and corrected release
fingerprints. Reuse the50-input public intake; two official mortality PDFs stay
explicitly prior-denied/not-retried under31's committed amendment. Require all25
data hashes, source fields, official URLs, current two-reader agreement and
full source marginal checks. Export actual full-release assay/weight marginals
and raw/native weight sums, adding inspectability that31 did not export.

Full-to-eligible masked design inventory must be preserved within each part:
train strata1–28/57PSU, test strata29–43/30PSU. One early stratum has3PSU;
other strata have2. Cluster=(stratum,PSU), never PSU alone. Require no source
participant or masked-cluster overlap between periods. No relabeling or
singleton rescue. Assign inner folds from eligible early clusters before
complete-case filtering; retain fixed cluster identities in all resamples.

Report actual events/non-events and event-bearing PSU counts in each part.
Require at least20 positive-weight events AND20non-events in each common part
before fits; otherwise export a blocked-feasibility diagnostic without model
metrics. This is a conservative computational/reporting threshold, NOT a power,
EPV, calibration-adequacy or clinical-validity certificate. Passing it does not
remove imprecision from a single later two-year cycle or few masked PSUs.

## Fixed models, contrasts and early-only tuning

Keep routine, routine_b2m, routine_cys, routine_cys_b2m exactly as31:
11 fixed terms age/age²,sex,current/former smoking,SBP,BMI,HbA1c,eGFRcrlow/high,
log2UACR; optional log2 reported Cys/B2M. Preserve31 cycle-specific creatinine
calibration, units/missing/LOD rules and first-three-positive SBP research mean.
No new predictor, interaction, cutoff, transform, clipping or imputation.

Primary: routine+B2M minus routine; primary metric late weighted8-year Brier.
Secondary: routine+Cys+B2M minus routine+Cys. AUC/calibration accompany both;
no B2M-vs-Cys winner contrast or change of question after positive results.

For each arm use3whole-PSU inner folds within EARLY TRAINING ONLY, salt
`forward32:inner:v1`; lambda1/10/100 minimizes pooled weighted inner-heldout
Brier, ties1e-12 prefer larger lambda. Weighted mean/population SD learned only
in each actual early training split; constants scale1. After tuning, fit fresh
preprocessing/coefficients on all early common cases, unpenalized intercept and
preserved ridge/Newton rules. Fit API must not accept late outcomes or weights.
Late feature changes may change predictions but cannot change early parameters.
Persist the four early-fit parameter artifacts as new research artifacts only.
Inner scores/optimization are training diagnostics, not new validation metrics.

## Point scores and two distinct uncertainty modes

Evaluate unchanged late predictions using native late weights: Brier,AUC,
observed/predicted risk,bias,O/E,offset calibration and intercept+slope,bins.
Calibration diagnostics never feed back into predictions. Report training and
late n/events, masked PSU/strata and coefficient/convergence/tuning diagnostics.

Exactly500 attempted paired draws, no redraw after failure. Separate streams:
late RNG seed20261103, early RNG seed20261104. Each iteration samples m_h−1
whole masked PSUs within each eligible stratum with factor m_h/(m_h−1), using
the native weights for that part. The SAME late draw is used in both modes and
all four arms; the early stream is independent. Preserve fixed early inner
folds; duplicate sampled clusters retain fold identity.

1. **fixed_training**: keep the single point early-fitted artifacts fixed and
   resample only late evaluation weights. This excludes training uncertainty.
2. **joint_refit (PRIMARY uncertainty)**: independently resample early weights,
   relearn all early preprocessing/tuning/full fits per arm, then evaluate on
   that iteration's common late draw. This reflects early fitting plus late
   sampling under the declared fixed design, not calendar/assay independence.

Keep separate results/CSVs/percentile2.5/97.5% bounds; never pool modes, select
the narrower favorable interval or call the fixed-training interval a refit CI.
Model/contrast failures and metric-local undefined calibration/AUC are retained;
an undefined secondary arm does not erase a defined primary contrast. Defined-
only bounds are not full nominal95% intervals. No training-uncertainty claim
from fixed mode, no assay/storage/selection-history uncertainty in either mode,
and no claim that this nonlinear refit bootstrap is an official NCHS variance.

## Assay drift, interpretation and delivery

Same public SSCARD assay platform is not proof of no drift. Baseline cycle,
storage duration and assay batch/order are not separately identifiable: public
release has no batch/run/individual-assay-date/storage-duration IDs. SSPRIS is
not a batch identifier. No invented batch correction or Cys standardization.
Report prespecified outcome-free mean/SD and missing/flag/pristine summaries
by temporal part; do not turn them into causal drift claims or subgroup tuning.

Executable tests cover native-period weights, split/source/design/support gates,
early-only preprocessing/tuning, invariance of models to late features/outcomes,
paired modes/RNG independence/PSU3, failures and exact aggregate exports.
Independent numerical/source audits, scientific CI and separate draft PRs.
Only aggregate parameters/metrics/source manifests/bootstrap draws, no source
records/participant IDs/prediction vectors. No accounts, agreements, restricted
data, contacts, merge, approvals or deployment.
