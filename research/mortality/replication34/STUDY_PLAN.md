# Iteration34: independent-replication readiness and precision, not model search

Written before executing the new planner or simulation cases. Preserve33B at
27ea33fc6cfd594c309b63ae59e695129511a9c8, neutral aggregate952a23344a6d16a453654de50bc4af437da605df9f8b2d120ced58fbbc9ec855.
All eight marginal deltaBrier/deltaAUC intervals included zero. No reinterpretation,
split search, new marker selection, model fitting or new outcome comparisons.
Laboratory reportv3/presentation, old models/results/main/UI are out of scope.
Joint integration and Fine–Gray remain deferred. No accounts, agreements,
restricted data, contacts, merge, approvals or deployment.

## Focused source qualification

Start from freeze30's named releases, current33A cause/assay/endpoint contracts
and programme-use history. The new question needs Cys+B2M, **not BTP**. Qualify
specific releases/visits, not an entire cohort name. Distinguish anonymous open
person-level releases, metadata/codebooks, aggregate publications, controlled
applications and provider prohibitions. An 'Open Study' catalogue label is not
permission to download participant data. Unknown evidence fails qualification;
it does not establish biological absence or global absence of a dataset.

For algorithm replication, documented interpretable clinical concentration
assays can differ from Siemens: do not impose the frozen30 CEP/BTP/brand gates.
NPX/RFU/DNAm surrogates are not absolute Cys/B2M concentrations for this question.
Numeric coefficient transport has a stronger mapping/artifact requirement.
No full33B final scorer exists: its40OOFfit records are not one frozen model.

Both underlying targets must map to cancer002/malignant neoplasms and broad
heart001, not cancer incidence, contributing flags, CHD alone or allCVD. Public
001/002/010 coarsening is sufficient to identify these two targets if010 is
known-other; all ten NCHS groups are not required here. Unknown early cause is
not known-other. Require exam origin,96months and known competing/vital outcomes.
Actual early censoring requires a separately precommitted IPCW extension; it
cannot be repaired by dropping censored people or inflating recruitment.

Preserve routine11 semantic definitions, same-visit requirements, appropriate
source-specific design/weights/families, interpretable units/QC/quantification,
joint domain and population support. Do not move an old urine/HbA1c measurement
to a later blood visit. Distinct recruitment and a programme history review are
required; no reidentification to prove zero person overlap. Previously used
NHANES/repeat participants cannot become untouched external data by a new file
name or assay/publication year. If no release qualifies, export precise
blockers/statuses/URLs and human next actions, without taking those actions.

## Precision scope fixed before calculation

For a paired fixed-rule Brier difference D=L_new-L_reference in[-1,1], assume
mu=E(D),variancev,actual evaluable paired n and a TOTAL loss-specific DEFF_D.
Conditional normal SE=sqrt(DEFF_D*v/n),halfwidth=z*SE,
n_precision=ceil(z^2*DEFF_D*v/h_target^2). Validate v<=1-mu^2. Use paired SD,
not marginal model-score SDs without covariance. Anchored SE needs reference
n andDEFF. Effect alone does not determine variance or required precisionN.
Never import33B effects/SEs as planner defaults or calculate observed power.

These are **conditional approximations**, not guaranteed sample sizes, clinical
acceptability, study power or adequacy for nested refitting. In primary
algorithm_replication mode always report algorithm_precision_unqualified:
training/tuning dependence is not captured by a fixed-loss formula. No default
training multiplier pretends to fix it. Source-coefficient transport is separate
and requires eight fully specified pre-external-outcome source scorers first.
The method proposal can be frozen now without pretending that scorer exists.

Illustrative numerical halfwidths/effects/incidences are technical assumptions,
not clinical thresholds or estimates from33B. LOW/CENTRAL/HIGH rows are scenario
ranges, not confidence intervals or probabilities of truth. Baseline incidence
is the **evaluable-sample event probability** for raw-count planning, not an
unqualified weighted national risk. Incidence alone cannot set Var(D).

Recruitment ceil(n_precision/(q_joint*q_known_conditional)) is expected-yield
accounting only. q_joint is availability of the entire qualified core+B2M+Cys
panel, not a product of marginal assay rates. q_known is conditional on that
availability. No implied MCAR correction or guaranteed enrollment yield.
q_known already includes endpoint-knownness losses; do not multiply another
(1-c) factor into that same yield. c is a separate methodological stop flag.
If before-horizon censoring c>0, flag IPCW_extension_required and disable current
uncensored-CIF eligibility regardless of a large recruitment number. Censoring
and linkage failure can be informative and change both estimand and variance.

Under explicitly independent Bernoulli events, give actual-count binomial
predictive quantiles at .025/.975 for n, not n/DEFF. Event DEFF_Y is separate
from loss DEFF_D. For unspecified clustered events, exact binomial ranges are
inapplicable; a scalar design effect is not a full event-count distribution.
No50-event or other illustrative count is labelled sufficient power. AUC and
calibration need separate precision targets/distribution assumptions.

## Synthetic numerical validation

Fixed finite prediction atoms define probabilitiesq,true risksr and paired
predictions beforeY. Compute exact incidence, paired mean/variance and simulate
only that hypothetical sampling law. No outcome-derived predictions or clipped
Gaussian loss. Equal-size clusters use a declared mixture: with probabilityrho
one common(Y,D) pair is repeated within the cluster, otherwise iid member pairs.
The exact marginal loss/event DEFF is1+(m-1)rho. Compare MonteCarlo moments with
the exact theory and report MonteCarlo uncertainty/coverage. No seed or tolerance
selection after outputs. This validates mathematical implementation under the
toy law, not assays, actual survey design, nested algorithm variance or clinical
performance. No synthetic result is external-validation evidence.

## Frozen method proposal and execution

REPLICATION_PROTOCOL.md freezes the two outcomes, four arms, contrasts and
methodological boundaries. It is a blocked proposal/template, **not a completed
cohort-specific preregistration**: a named qualified release, exact mappings,
appropriate clustering/censoring strategy and investigator-selected precision
targets must be fixed before new outcome inspection. None are inferred from a
favourable old effect. Primary mode is algorithm replication; locked transport
is separately gated. Do not fit or score any new cohort in34.

Implement a pure qualification gate, planner and aggregate CLI with explicit
input schemas, source/review hashes and preserved103prior scientific modules.
Freeze SCENARIOS.json/SIMULATION_CASES.json before their production execution.
Tests cover algebra/bounds, anchoredSE, DEFF double counting, jointavailability,
exact-binomial boundaries, censoring/unknown fail-closed rules, algorithm vs
locked-transport claims, deterministic bounded simulation and aggregate exports.
Scientific CI runs the scenario package without participant-source downloads.
Independently audit arithmetic/intervals/source classifications/code/provenance;
export only assumptions, metadata and aggregate simulation summaries. Separate
draft PR/evidence, no changes to33B or any laboratory artifact.
