# Iteration33A: source and event gate before separate-cause marker comparisons

Written before new cause-specific counts, weighted rankings, fitting or comparative
metrics. Iteration32 is saved at fe4d20672dc761c8a6533b3dcaf5a6a0340251d8
(draft theoriesPR24; evidence1551893298feff1bd7ce41e97d06c1faf6af60a1, PR125).
This is the next stage of the existing research task, not a scheduled task.

## Question and boundary

Determine which of up to five leading **defined** public death-cause groups can
support a separately specified marker comparison in this cohort. Rank by weighted
eight-year mortality burden, accompanied by actual event counts, never by AUC.
The ranking describes this received NHANES domain, not Russia, the world, or an
unqualified US population. Heterogeneous UCOD010 is reported for accounting but
is not a single clinically directed cause. Do not force five eligible groups.

This stage exports source/measurement/time/weight/design/cause counts and support
diagnostics only. It does not fit models, tune parameters, compute predictive
metrics or bootstrap model uncertainty. A separate immutable comparison protocol
must precede any subsequent comparative fits. The gate is not itself a marker
benefit result or an independent validation. Joint integration and Fine–Gray
remain deferred; neither is implemented or run here. Historical five-cause
piecewise exponential ridge/CIF and all-cause models remain unchanged.

## Source, domain, time and measurement contract

Reuse the exact50-input official-source intake and25 historical data SHA pins of
b2m31, its two-reader mortality check and corrected July2022 SSCARD_A release.
No fresh attempts at the two previously denied mortality PDFs, ZIP downloads,
alternate mirrors, accounts, agreements or restricted data. Record both denied
PDFs as prior-denied/not-retried, separate from50 successfully retrieved inputs.
Existing retained endpoint_contract11 and official-source review define public
initial underlying UCOD_LEADING groups001..010, not contributing-cause flags.
The code/time release is2019 linked mortality; official public cause/time
perturbation remains a limitation. Do not infer individual restricted causes.

Baseline1999–2000,2001–2002,2003–2004; original MEC exam is time zero, not the
2018–2020 serum assay date. Inherit b2m31 eligibility and the same primary common
flag00 domain: age40–75, examined, linked eligible, known mortality status and
nonnegative finite exam followup, positive corrected pooled subsample weight,
all11 routine terms and both positive finite reported proteins with both flags0.
No new exclusions based on a cause or a protein's predictive association.
Known prior common n6044/events8=586 are domain consistency checks, not new cause
results. All cycles/participants have prior programme use; no new holdout claim.

Endpoint horizon is exactly96months, inclusive; death at0 or96months counts.
Every eligible survivor must be observed through96months. Deaths after96months
are known noncases by this horizon. Require known public cause001..010 for
every common-domain death at/before96months before allowing any later cause
comparison. Unknown early causes block, are reported separately, and are never
silently treated as negative for all causes. Missing cause after96months does
not make an eight-year label unknown but is reported. A cause recorded for an
alive participant, fractional/undeclared cause, unsupported cycle, invalid
status/time, changed source or design is a structural error and stops the gate.
Never substitute a coarsened label or a subtype (e.g. allCVD for heart001).

All three baseline cycles expose ten public groups; the2015–2018 reduction to
001/002/010 concerns baseline survey years, not death years. Refuse later cycles
instead of assuming absent groups have zero deaths.

Pooled six-year estimand, **not32's temporal estimands**: (2/3)WTSSCB4Y for
1999–2002, (1/3)WTSSCB2Y for2003–2004. Verify this formula against the retained
native fields. Common complete-case/measurable coverage, missingness, source
flags, assay calibration/storage/QC/LOQ uncertainty remain; weights do not undo
selection. No new IFCC correction, III calibration, imputation, clipping, tail
threshold or assay harmonisation. Cys/B2M are reported serum mg/L, log2 in any
later fixed comparison; this gate does not estimate their effects.

## Counts, rankings and predeclared support

Export all ten groups plus unknown cause: overall/percycle n, eight-year events,
positive-weight event rows, event-bearing masked PSUs/strata, weighted event total,
weighted eight-year incidence, share of all eight-year deaths, and source-cause
counts at all observed followup times. Also report squared event-weight sums,
event-weight Kish effective n and maximum single-event weight fraction; these
are descriptive concentration diagnostics, not additional selection thresholds.
Common domain denominator is identical
for all causes. Account for all586 prior all-cause events and their weight sum
using the disjoint ten groups plus unknown. Keep010 visible, out of clinical
ranking. Weighted incidence is a burden estimate, not a fitted prediction.

Rank001..009 by decreasing common weighted event total, then decreasing raw
event count, then code ascending to resolve exact ties. The first five are
priority candidates regardless of their support; do not replace an unsupported
top-five cause with a lower-ranked convenient one. Report all nine ranks and
every failure. At most five, possibly fewer or zero, proceed to a later protocol.

Retain the outcome-independent b2m31 whole-masked-PSU fold maps: five outer
folds, salt b2m31:outer; three inner folds within each outer training, salt
b2m31:inner:{outer}. Construct outer map on the full eligible design before
common filtering, preserve all87 masked PSUs including stratum1's third unit,
and require the same eligible/common inventories. Report actual label support
in every outer and inner training/evaluation cell, with no search for a better
partition. Masked PSU IDs are public resampling units, not participant IDs.

Predeclared pragmatic minimum for a top-five cause to proceed:

- at least50 positive-weight event rows and50 positive-weight noncase rows;
- events present in at least20 masked PSUs;
- every outer evaluation has at least5 event and5 noncase rows;
- every outer training has at least20 event and20 noncase rows;
- every inner training has at least10 event and10 noncase rows;
- every inner evaluation has at least3 event and3 noncase rows.

These conservative reporting/computational thresholds are not universal power,
EPV, identifiability, calibration or clinical-readiness certificates. They do not
guarantee defined bootstrap or calibration diagnostics. Keep rare causes and
their precise failed cells visible instead of combining biologically distinct
causes or weakening thresholds after counts. Passing means feasible for an
exploratory comparison only. Count-gate selection based on outcomes is explicit;
it does not establish confirmatory inference for the selected causes.

## Contract for a later separate comparison

If source/cause/no-censor gates pass, a fixed96month binary target
I(death<=96 and cause=k) can estimate that cause's cumulative incidence directly.
An earlier competing death is a known noncase for this target, not a censored
observation and not an assertion that the person is alive. A later death and a
survivor observed through96months are also noncases at this horizon. Do not use
this construction when censoring or early unknown causes make labels unknown.
This is not a cause-specific hazard, Cox model or Fine–Gray model. Independently
fitted probabilities must not be summed or advertised as a coherent joint set.

A later protocol will specify routine±B2M±Cys comparisons on the identical
participants/folds/weights for each supported cause, training-only preprocessing
and nested tuning, Brier/discrimination/calibration and paired full-refit
uncertainty with metric-local undefined outcomes and multiplicity disclosed.
No marker is assigned exclusively to a single cause. No benefit, optimal panel,
causal effect or clinical validity follows from this source/count gate.

## Reproducibility and exports

Commit this plan, cause definitions, support rules and SHA lock before executing
the gate. Add a new wrapper importer and tested pure gate without modifying
earlier modules/results. Pin all95 prior scientific modules, old source/assay
and endpoint contracts and frozen artifact. Run synthetic source/endpoint/
weight/competing-outcome/rare-event/fold tests on Python3.11/3.13 and the actual
official-data gate in scientific CI. Export only aggregate JSON, exact source
manifest, code hashes and event/fold/PSU summaries; never participant IDs,
individual labels, predictions or raw source files. Independent audit should
recompute count-accounting, weights/ranks/support rules from exported aggregates
and compare stdout block hashes; it must disclose no independent raw reanalysis.
Use a separate draft PR, verify remote commit and CI, never merge or deploy.
