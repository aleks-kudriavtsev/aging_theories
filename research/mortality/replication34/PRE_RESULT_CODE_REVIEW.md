# Rep34 final pre-production code review —2026-10-29

Disposition: **NO remaining scientific/code blocker found for the declared
metadata/scenario/toy-law stage.** This is not approval of a new cohort or of
refitted-algorithm precision, clinical validity, data access or model deployment.

No repository edit, GET, new outcome/participant read, production scenario
calculation or production20k Monte Carlo case was performed in this review.
Only code, committed method/scenario inputs, tests and metadata were read.
The metadata-only check_study function was run independently:22 protocol
metadata SHA pins,103 historical scientific Python SHA pins, SOURCE_LOCK and
source review fingerprints all passed without planning/simulation calls.

## Reviewed fixed scope

- Protocol4df636ea18ca8bb29b05788c37800438566426ec predates the12 stipulated
  scenarios and4 production toy simulations. Fields/rows/seeds remain fixed.
  Runner passes each typed scenario directly to precision.plan; it does not
  parse33B deltas, SEs or bootstrap results. Old33B files are read only as bytes
  for preservation hashes, not as effect/variance defaults.
- Source gates are pure metadata/citation qualification, not person-level data
  loading. They distinguish algorithm replication from locked-source transport;
  processing/download/outcome activation remains false. Transport requires a
  separately frozen source scorer and assay mapping; algorithm-level precision
  remains unqualified. No model/outcome pipeline is imported or called.
- precision.plan accepts only its explicit schema, one total DEFF basis and one
  paired variance source. Anchored SE needs N_ref/DEFF_ref. Loss DEFF and event
  DEFF are separated; exact-binomial predictive quantiles use RAW n only under
  stated event independence, never n/DEFF. Recruitment is expectation-only
  conditional availability accounting, not assurance or selection correction.
- All plans retain full_CIF_validation_eligible=False, algorithm adequacy
  unqualified, no AUC/calibration/power/clinical-adequacy certification. Any
  positive pre-horizon censoring flags an IPCW extension, not silent deletion
  or negative labelling. No independent cause probabilities are summed.

## Numerical boundary corrections independently re-read

The positive-variance coverage rule has NO added absolute slack. The earlier
1e-14 problem, which would turn a valid tiny-loss exact87.5% known-z law into
100% coverage, is removed. Only exact-zero-variance laws use the explicitly
scoped16machine-epsilon*loss-magnitude constancy check of participant losses
and their mean; it cannot widen positive-variance intervals.

Moment feasibility now uses stable bound (1-|mu|)*(1+|mu|), permits only8eps
RELATIVE arithmetic roundoff with abs_tol0, and rejects positive-scale variance
underflow. Thus mean±1 cannot admit positive SD/anchored-SE variance and clip it
tozero. Near-endpoint and SD/SE unit fixtures cover this. Neither correction
changes the committed input DGP, seed, scenario, target or observed model result.
PRE_RESULT_IMPLEMENTATION_REVIEW.md records both pre-production findings.

Prediction-atom moments respect X→fixed paired predictions→Y, with exact total
variance components. Shared-pair cluster mixture has declared equal-size-cluster
DEFF1+(m-1)rho; it validates only that toy law. Centered fourth moment is exported
so the approximate Monte Carlo variance-SE can be independently recomputed.

Runner checks meanD,varianceDbar and mean event count using the declared6MCSE
band (three checks per each of four cases); zero MCSE uses exact equality.
CI coverage is reported with MC uncertainty and is NEVER required/tuned to95%.
Rare/skewed and degenerate laws can legitimately show different normal-interval
coverage. Successful toy checks do not validate a refitted algorithm, external
cohort or actual calibration/discrimination precision.

## Export, workflow and tests

Exactly five aggregate files: PRECISION_RESULTS.json,SIMULATION_VALIDATION.json,
READINESS_RESULTS.json,SOURCE_MANIFEST.json,RESULTS_SUMMARY.json. Source manifest
bytes are copied exactly, four payload hashes enter the summary, and every file
has an exact REP34_FILE_BEGIN/END stdout block. Complete serialization and a
second metadata check precede any output-directory creation; failure leaves no
partial success export. Toy disagreement is retained in diagnostics and fails
the portable CLI rather than being removed or rerun with a favourable seed.

Workflow installs only NumPy/SciPy and runs the pure/mocked41-test suite on
Python3.11/3.13; production scenario job is explicitly separate on PR/Python3.12.
Artifact upload contains only the aggregate output directory. No source-data
retriever, account, agreement, restricted-provider access, merge or deploy step.
Runner tests replace BOTH precision.plan and simulate_case with mocks; their
12/4 call checks do not execute production inputs. The reported41 passing unit
checks are consistent with the reviewed tests; this review did not rerun them
or duplicate production calculations.

Reviewed current SHA256:
- precision.py: a40d8a79ff4428c674795f32bb4aa63998d8912546aaa473fc6e537a99486bcc
- gate.py: 7ac8f6f5544254c8b828387805470128f5bb03e7500f4f8023b2756bc5c8a9b7
- analyze.py: 3dd2e087555770c5bba729f4ee04e2b4bbc47aad40d3d98466ea40458b0f8d38
- PRE_RESULT_IMPLEMENTATION_REVIEW.md: 549f5876134c08044265581688e7386b679421e7a8d641e674a7857bc01a43fe

Source/runner scope was independently reviewed in parallel without production
execution; final notes are available separately from the same review thread.

Parallel independent runner/source review completed with no blocker:
/tmp/rep34_runner_source_review.md. It confirms the same analyze.py SHA and
metadata/source pins, no33B effect parsing, five exports/markers and no
coverage-based pass criterion. One conservative generic limitation is retained:
with exact-zero MCSE, strict equality can reject tiny roundoff for a nonzero
constant toy law; that would fail a check rather than manufacture success or
95% coverage. It does not justify weakening checks or changing frozen cases.
