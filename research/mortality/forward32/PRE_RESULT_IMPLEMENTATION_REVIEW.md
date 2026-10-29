# Iteration32: independent pre-result implementation review

Date: 2026-10-29. Reviewer: delegated agent `/root/assay28`.
Scope: committed protocol377778ca0ddbf48af5e0dc1735b88c24c6b8a11c, LOCK,
native-weight rules, new model/cohort/retrieve modules, and inherited numerical/
source qualification code. No participant source records, temporal event counts,
new comparative metrics or job logs were read. No network requests or data/model
calculations were performed. Repository code was not edited by this reviewer.

Final status: model/source and final orchestration/bootstrap review completed
before any new temporal results were read or calculated. No remaining critical
scientific implementation blocker was found in the reviewed code. The initially
identified fold discrepancy was fixed before metrics, as documented below.
Scientific CI and synthetic tests are separate validation by the implementation
team; this review does not claim their eventual result. This is a scientific
code review, not PR approval, assay certification or clinical readiness.

## Protocol and provenance checks

The committed ANALYSIS_PLAN SHA256 is
`9eb7ca3e3faf076bb34c1c1e42bd6cb4f251ad7aac53165d07f57bae578a0cd4`, matching
LOCK and the current plan. Calling only `retrieve.check_study()` initially
verified13 metadata fingerprints and90 historical scientific Python fingerprints.
After the pre-result implementation notes were committed at
`0a277179fbf3104f51086ad8899132e13d96071e`, I reran the same check: all14 metadata
fingerprints and all90 historical Python fingerprints pass. This reads
local source/metadata bytes for hashing only and does not load participant tables.
The source API declares48static files plus catalogue and resolved SSCARD=50inputs.
The previously denied two mortalityPDFs remain explicitly `denied_not_retried`;
there was no new PDF request or attempted workaround in this review.

The design fixes early biological baseline1999–2002 for learning, late2003–2004
for evaluation, the same11 engineered routine terms, four inherited31 arms,
primary B2Mminus routine contrast, secondary B2Mafter Cys contrast, and8-year
binary all-cause mortality. Biological baseline is exam, not retrospective assay
year2018–2020. Persons and the later calendar period were already used in this
programme: this is exploratory temporal stress evaluation, not an untouched or
independent external validation and not transport of frozen30 models.

## Source, measurement, weights and support

The temporal adapter explicitly selects raw corrected WTSSCB4Y for early and
WTSSCB2Y for late, keeping W6 separately for provenance. It does not require late
WTSSCB4Y, which is structurally missing. It does not normalize the early cycles
separately. The raw selected native weight is used for scoring; likelihood
normalization to actual training n is a separate declared ridge convention.

The effective target is each period's received complete-case/measurable-protein
and linked domain. Native survey weights do not automatically correct additional
linkage, item-missingness or both-flag0 selection. An unqualified national claim
would be unjustified. The code exports native/rawW4/rawW2/W6 weight support,
Kish n and complete-domain weighted coverage separately by part and cycle.

The unmodified pinned31 loader supplies exact clinical fields/units, corrected
creatinine mapping, primary both-flags0 common mask, positive finite reported
Cys/B2M, deterministic XPT zero decoding, historical25data-byte pins, full-release
marginal guards, exact keyed source joins,2019public LMF and two mortality-reader
agreement. The new adapter preserves common-mask identities, checks full eligible
masked inventory separately for strata1–28 and29–43, no shared participant or
masked PSU identity, no singletons and native positive finite weights. Public
masked inventory, rather than PSUcode alone, is used; original early PSU3 is not
silently forced into a two-PSU convention.

Reported protein values remain assay/platform-specific prognostic signals.
Detectionflag0 is not per-specimen QC or established LOQ validity. No unpublished
batch identity, cystatin traceability correction or storage-time correction is
invented. Baseline calendar, storage duration and assay batch/order are not
identifiable separately in the released variables. SSPRIS is not a batch ID.
Neither bootstrap captures assay/storage/selection-history uncertainty.

The predeclared gate counts actual positive-weight event and nonevent rows in
BOTH common parts (minimum20each), with death after96months a nonevent. It checks
uncensoring before model fits. A failing gate must produce blocked-feasibility
aggregates without fitting or invented metrics. This minimum is a pragmatic
computational/reporting rule; passing is not a power, EPV, calibration-identifiability
or clinical-validity certificate. Event-bearing masked PSU counts still matter.

## Pure model boundary and early tuning

`fit_forward` and `fit_training_arm` do not accept late outcomes or weights.
Late features are validated and scored but do not enter preprocessing, penalty
selection or coefficient fitting. Every arm's3inner training splits independently
learn weighted means/populationSD; the full early fit learns them afresh after
pooled early inner-heldout weighted Brier selects lambda1/10/100, tolerance1e-12
preferring larger lambda on ties. The intercept is unpenalized. Full-fit artifact
contains coefficients, preprocessing, convergence/inner support aggregates, and no
participant identities or probability vectors. `predict_forward` only replays
fixed early parameters; there is no late intercept or recalibration update.

One pre-result discrepancy was reported to parent: `_validated_split` initially
called `psu_partition` on common early clusters, although the plan specifies folds
assigned from the eligible early cluster inventory before complete-case selection.
The inherited partition sorts cluster hashes and assigns rank%3, so disappearance
of a whole common-case PSU can alter other folds. Required resolution before fits:
pass/use the eligible-derived fixed mapping, or explicitly prove/gate complete
eligible/common cluster equality. Bootstrap zeros must retain cluster/fold identity.
Resolved before metrics: the temporal source gate now requires common masked
inventory to equal eligible inventory in EACH partition and exports the actual
early eligible-derived threefold mapping using the locked salt. Equal cluster
sets imply the identical hash-sort/rank partition, while zero-weight retained
rows preserve identities in resamples. Any lost common PSU produces an explicit
blocked-feasibility reason before fits, rather than silently changing folds.

## Final orchestration and paired uncertainty review

The reviewed `forward32/analyze.py` keeps source/schema/design/support checks
before point fitting. The blocked-feasibility branch exports manifest/source and
actual reasons with0model fits and0draws, then returns. It cannot produce model
metrics from a blocked source domain. Native period weights must exactly equal
the appropriate retained raw W4/W2 values before arrays are prepared.

For an allowed gate, early and late arrays are extracted on their respective
fixed common masks. The common/eligible PSU-set equivalence is rechecked, and
the recomputed model fold mapping must match the source audit's eligible-derived
mapping. `_check_inner_mapping` and `_fit_arms` pass no later outcomes/weights to
model fitting. Point arms are fitted separately through `fit_training_arm` after
structural checks, so an individual fitting failure is recorded locally rather
than erasing another valid primary pair. This also resolves the potential
all-or-nothing behaviour of the pure convenience `fit_forward` helper.

`bootstrap_modes` constructs separate RNG objects with locked late seed20261103
and train seed20261104. On every iteration it calls the inherited rescaled
whole-PSU sampler once on the full eligible late frame and once on the full
eligible early frame BEFORE selecting common rows. That sampler draws m_h-1
units with replacement and multiplier count*m_h/(m_h-1), supporting original
three-PSU strata. Both modes and all four arms use the SAME late weight array.
The early array is used only in joint refitting; the fixed mode reuses the
original point predictions and never fits. Keeping zero-weight rows prevents
cluster/fold reassignment in resamples. Streams are not reseeded, redrawn or
advanced conditionally after fitting failure. Exactly500attempted draws per
mode are requested from the LOCK-checked configuration.

Each joint arm independently relearns inner preprocessing, chooses lambda on
early inner-heldout Brier and relearns final preprocessing/coefficients. Arm
failures are retained; contrasts and metric-local AUC/calibration undefinedness
are handled independently. A missing secondary model cannot erase a defined
primary contrast. Selection counts correspond to successful training fits,
while separate evaluation/contrast statuses describe score availability.

Native late weights are passed unchanged into scoring in both modes. The
inherited score functions do not mutate predictions or weights; their clipping
is confined to diagnostic calibration logits. Late offset/intercept/slope
estimation never updates scored probabilities. Weighted quintile calibration
bins are aggregate diagnostics, not a late fit/calibration step for the model.

Mode-specific CSVs retain all attempted rows and statuses; paired contrast CSV
retains each mode/replicate identity. Intervals verify total-row counts and
report valid/undefined counts per metric, with scope relabelled to the actual
fixed-training or joint-refit procedure. joint_refit is explicitly the PRIMARY
uncertainty; the modes are not pooled. Defined-only bounds are not guaranteed
full nominal95% intervals. A fixed-mode interval omits learning uncertainty.
Neither is an official NCHS nonlinear-model variance estimator or a measure of
unavailable batch/storage/marker-selection-history uncertainty. Distinct masked
period labels and independent RNG streams do not prove that true geographic
sampling units were independent across periods.

Exports contain source metadata, source/design/feature summaries, model
parameters/inner support, aggregate late scores/calibration bins, bootstrap
aggregate rows and hashes. Participant records, IDs, predictions and raw
source measurements are not exported. Partial point artifacts carry explicit
per-arm fit statuses; downstream scoring helpers correctly demand a complete
four-arm artifact rather than silently treating absent coefficients as valid.

Reviewed-code SHA256 at final read (all files new to iteration32):

- model.py: `d1389c92c0276698dfbae0a27692a233ce1c0bd0afb987b02384db77ae0a511f`
- cohort.py: `67e45f7a4e3564842ef1ecdc3a6011b1c70f64dac36563e3f62e533cc5d050fb`
- retrieve.py: `262effbe5e544929227f55560ac1eb1dd5930ea600d92a4d158f79e854b7caa5`
- analyze.py: `f956eb363c9172fecb9e2f865051623596290fb3f8d8bbecfc83ea948bfda487`

The source/model/scoring code was reviewed without running an empirical intake
or fitting models. Final tests/scientific CI remain the implementation team's
pending check at this review boundary. No new denied, pending or uncertain
network/approval action was attempted by this reviewer. Historical mortalityPDF
403denials remain explicitly denied_not_retried, not new unfinished requests.

## Pre-result interpretation clarification: label-availability time

The split is retrospective by BASELINE calendar, not a historical as-of2003
simulation. Eight-year labels of early1999–2002 participants become observable
in2007–2010, after late2003–2004 baselines; the current endpoint uses2019public
LMF. Likewise both early and late baseline sera were assayed in2018–2020.
Excluding later-participant outcomes from fit is scientifically useful, but does
not imply that these trained models or complete early labels were available when
the late participants were examined. Do not describe this as a historically
prospective2003 deployment validation. No horizon, algorithm or contrast change
is needed: retain the declared exploratory retrospective baseline-calendar target.
