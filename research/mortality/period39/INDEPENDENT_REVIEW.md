# Period39 independent method review — pre-output boundary

2026-10-30. Worktree `/workspace/aging_theories_period39_local`; frozen protocol
`c0a570e0ac40de94da8a05b53f5e9d80aaf7aacf`. This initial review covers the protocol,
study plan, source manifest and cached extractor/docs. `reference.py`/runner were
not yet present when this note was started. No actual-source probability product,
new source GET, original RDA, participant data, fit, outcome comparison, producer
CLI, repository change or external write was performed. A later code review and
independent numerical audit remain necessary after the corresponding freeze.

## Initial decision

No critical methodological blocker found for a **conditional local-derived
population period reference** under the explicitly declared units and schedule.
It must not become an authenticated official WPPqx release, an individually
calibrated forecast, a2026–2036 calendar projection, a disease-free/stage-specific
baseline or biomarker model. Source raw-RDA scaling and original upstream bytes
are not independently reauthenticated; this is a material qualification, not a
statement that original source validity has been established by local hashes.

## Arithmetic and interpretation

1. For exact birthday age a and h whole years, use annual intervals with attained
   ages a,...,a+h−1. Baseline80, horizon10 ends at birthday90 and needs80–89,
   not an age90 rate. A period schedule fixes the same calendar2023 age-specific
   rates at every future attained age; no future-calendar rate path or cohort
   trend is identified. Alive-at-exact-birthday is a synthetic period condition,
   not a measured state or health eligibility of an individual.
2. Primary CF uses `Q=1-exp(-sum mx)` and annual survival `exp(-mx)`. This is the
   declared constant-force approximation. Central death rate mx is not itself
   a one-year probability and must not be multiplied as if it were qx. Annual
   rate units must be deaths/person-year, not per1000, percent or cumulative
   hazard without a time interval.
3. Secondary UDD with ax=.5 uses `q=mx/(1+.5*mx)` and survival
   `(1-.5*mx)/(1+.5*mx)`. It is a different lifetime-within-interval assumption,
   not a data-selected alternative. Its valid probability domain is0≤mx≤2;
   values>2 produce q>1 and must fail, including any rounding upper endpoint
   outside that domain. No clipping/renormalization should repair such inputs.
4. Both constructions satisfy chain survival exactly in mathematical arithmetic:
   `Q(a,10)=Q(a,5)+(1-Q(a,5))*Q(a+5,5)` when the **same** attained-age schedule
   and method is used in each term. Using one age repeatedly, different calendar
   years, a CF/UDD mixture or independent population samples would not be this
   identity. Ten-year risk cannot fall below five-year risk under nonnegative
   hazards and valid annual probabilities; age-order monotonicity need not hold
   and should remain a diagnostic, without forced smoothing.
5. There are41ages×2sexes×2horizons=164 reference cells. There are36start ages
   40–75×2sexes=72 chaining rows; if both CF and UDD are shown, the report should
   clarify that each row checks each method rather than imply only72 scalar
   method checks. Starting76–80 cannot chain through an81–85 baseline entry in
   the predeclared164-cell output, despite the underlying age schedule coverage.
6. F/M are source categories, not an automatically chosen person category. Sexes
   cannot be averaged without population weights; a local schedule cannot
   establish representativeness for other sex/category definitions.

## Six-decimal quantization and stored qx

The pinned cached extractor uses `round(mx,6)` and independently computes
`round(1-exp(-original_mx),6)` **before** exporting rounded mx. Therefore stored
qx is author-derived, not an official observed WPPqx; exact equality with CF
recomputed from rounded mx is not expected. A1e−6 discrepancy tolerance is
consistent with the two half-unit decimal-rounding errors, since
`|d(1-exp(-m))/dm|=exp(-m)≤1` for m≥0. This does not verify upstream scaling or
make a discrepancy inside tolerance evidence of demographic calibration.

Under the explicit round-to-nearest6-decimal assumption, each underlying rate
lies in the closed interval `[max(0,m−.5e−6),m+.5e−6]`. Ties-to-even or endpoint
attainability does not invalidate conservative closed bounds. Propagate all
lower rates and all upper rates separately, using the same selected method.
These are arithmetic sensitivity envelopes, not sampling confidence intervals,
forecast intervals, or independent errors sampled from a probability model.
No independence distribution is required for a conservative box envelope;
calling the intervals “independent source_mx” must not imply such a distribution.

CF is coordinate-wise monotone with derivative `exp(-sum m)≤1`. UDD q has
`dq/dm=(1+.5m)^(-2)≤1`, and the derivative of the multi-year death probability
with respect to any annual q is a product of other annual survivals in[0,1].
Thus for valid UDD rates and CF, total envelope width is at most the sum of rate
widths, ≤h×1e−6. The protocol's small numerical allowance is a computational
guard, not epidemiological uncertainty or a statement that the two methods
span all realistic model error. CF-vs-UDD differences must be labelled a
predeclared assumption sensitivity, not a CI.

## Population mixture and double counting

The source schedule includes people with disease, prior events, treatment,
comorbidity and heterogeneous markers. Population all-cause mortality is already
a mixture. It cannot be labelled the intensity of healthy people, used as an
untreated counterfactual, or added to postdiagnosis all-cause death/CIFs. Such
addition generally double-counts deaths and cannot be fixed by clipping a sum
to1. Calibration of a conditional disease-state model to this anchor would need
representative state/risk-factor distributions, compatible attained-age/calendar
conditioning and independently checked outcomes. None is supplied here.

The strict adapter should require explicit exact-birthday integer age40–80,
RUS/2023, sourceF/M, horizon5/10, the period-schedule assumption, declared units
and `population_reference_only` purpose. Unknown fields, biomarker inputs,
disease-free/stage requests and unsupported sex/country/year/horizon should fail.
There is no permission to silently ignore a person's clinical information while
returning this reference as their risk.

## Source checks completed before parsing

All four current local source files and their Git blobs at
`31fff9249237a33d9a10ac107d6710a04bfec517` match the protocol hashes, including
CSV SHA`59104bce947c3ab1132893e71decf946af3725c151fa5ff6f64bef440d4b54a2`.
Machine-readable evidence: `/tmp/period39_preparse_source_hashes.json`.
No CSV rates or probability products were parsed/computed for this check.
These hashes prove local derivative identity, not original RDA content, source
licence authentication, numeric scaling or a latest official release claim.
2023 is the latest year in the pinned local estimate CSV only.

## Pending code/audit checks

After code is available, review exact schema and decimal-resolution checks,
nonnegative/finite rates, valid UDD domain including envelope endpoints,
qx diagnostic scope, explicit clock/purpose/unit gates, no extrapolation and
numerically stable CF/product formulas. After code freeze, an independent
stdlib verifier should check Git/current code and four source pins, all164
85-digit Decimal references,72 chaining rows for every advertised method,
rounding envelopes, horizon coherence, source/target/period metadata and
compact output hashes. It should reject3–5 deliberately corrupted aggregate
fixtures, import no producer code and run no producer CLI or new source fetch.
SHFM review remains separate, source-only and blocked until a complete qualified
numerical contract is fixed.

No action was denied, pending or uncertain in this review. No source-access
restriction was retried. The incomplete code review and numerical audit are
normal stage boundaries, not an external-access or authorization blocker.

## Pre-freeze code review — initial findings

`reference.py`, `analyze.py` and `plot.py` are now available and have been read.
The production path hashes all4source files before parsing; validates100 keys;
uses attained ages a through a+h−1; preserves immutable Decimal maps; performs
both methods, stored-qx products and rounding envelopes at85 digits; and exports
explicit period/mixed-population/source-unit limitations. UDD checks the upper
rounding rate≤2 and rejects saved rate2 because its envelope would exceed the
valid domain. The runner exports72 chaining age/sex rows with144 method checks,
164 horizon-coherence checks and clear non-CI/source-scaling flags. These
arithmetic/domain components agree with the frozen protocol.

Two adapter/presentation findings were sent to root **before code freeze and
before any actual164-cell probability calculation**:

- **F1, medium public-API correctness:** `score_reference` currently accepts a
  `ReferenceBundle(source_files_verified=False)` returned by the pure parser,
  yet returns a named `RUS_WPP2023_period_reference` result. Warning metadata is
  truthful but does not meet the strict source-specific adapter contract. The
  production runner is unaffected because it uses pinned `load_bundle`.
  Recommendation: reject an unverified bundle in this public adapter, preferably
  also require the exact source commit/hash inventory/unit-evidence identity.
  `calculate_reference` already provides the appropriately unnamed pure-math
  fixture API. Synthetic positive adapter tests may explicitly mock the
  verified loader; they must not imply new evidence for the real source.
- **F2, low presentation/source-label correctness:** `plot.series` initially
  checks164age/sex/horizon keys and a top-level scope label but ignores country,
  source year, period assumption and purpose in the request. Thus an otherwise
  similarly structured USA/2024 payload could be plotted with a Russia2023
  title. Current producer output is unaffected. Recommendation: require fixed
  top-level and per-row context and exact integer age/horizon types, without
  calculating any new probabilities or broadening the model.

No other critical arithmetic issue was found in the inspected code. This note
records open pre-freeze findings, not a claim that root has already fixed them.
No actual-source probability product or full producer scenario CLI has been
executed by this reviewer. Parent-reported25unit completions have not yet been
independently rerun by this reviewer.

Reviewed initial code SHA256:
- `reference.py`: `a1ab2c473d7ee20d273555a00bff577d0474191dfc4896ae219398b87cd8c6ce`
- `analyze.py`: `e23b3cb61ab92495e75c8425c8cbf49cf7da34ccc6e953675f7c065c8944157f`
- `plot.py`: `51edb85829b20bf4e8ea91817eefb67011e9bdb49414b8c457773a7b5801e686`

## Final pre-production follow-up — F1/F2 resolved

Root addressed both findings before computing the actual164-cell source grid.
Final `score_reference` rejects `source_files_verified=False` before arithmetic,
then requires source commit, exact4-file hash inventory and unit-evidence
identity to match the immutable manifest. Pure mathematics remains available
through `calculate_reference`; positive named-adapter unit tests use an explicitly
mocked manufactured loader, not a claim of authentic source rates.

Final `plot.series` requires the fixed Russia2023 exact-birthday period/purpose
context, full request schema, exact integer coordinates and all164keys. The
country/year relabelling gap is closed. The source-specific API is a fail-closed
research contract, not a security boundary against deliberate Python-object
forgery or proof of original-RDA scaling.

Both fixes have been read and independently tested locally with manufactured
fixtures only:23 core tests and4 analyzer/presentation tests, **27/27 PASS**.
No actual source probability grid or producer CLI was run by this reviewer.
There are no remaining critical scientific/code findings for the stated
conditional local-derived period scope. The later actual aggregate audit
still needs the frozen code identity, output bytes and logs from root.

Final reviewed code SHA256:
- `reference.py`: `a1ab2c473d7ee20d273555a00bff577d0474191dfc4896ae219398b87cd8c6ce`
- `analyze.py`: `e23b3cb61ab92495e75c8425c8cbf49cf7da34ccc6e953675f7c065c8944157f`
- `plot.py`: `51edb85829b20bf4e8ea91817eefb67011e9bdb49414b8c457773a7b5801e686`

## Final independent aggregate audit — completed

Final scientific freeze is `151fe91b1540b91926baf3b3cc782cc7515ae55b`.
The subsequent separate stdlib tool `/tmp/period39_independent_verify.py` imports
no producer module and performs no producer replay, source GET, participant
processing or model fit. Its compact JSON/Markdown audit reports
**10295/10295 PASS**, both against the original nested output folder and against
the preserved flat repository results. It verifies14current/Git code pins,
3unchanged protocol files,4current/Git source pins,3source-review hashes,
100target input keys, all164reference cells/CSV rows,20selected rows,
72chaining rows×2methods,164horizon checks, source-age80→89 support and all
scope/source/clinical flags.

The independent85-digit oracle uses a product of annual exponential survivals
for CF and exact rational survival products for UDD and stored-qx diagnostics.
Every probability, rounding endpoint and diagnostic is checked, rather than
only a selected table or maximum. Both chart renderings have identical SVG/PNG
bytes, matching their clean provenance log and source-result SHA; no plot or
producer calculation was repeated. Chart-byte authentication is not an
independent rendering/epidemiological validity check. Root reports its own
visual inspection separately.

Five intentionally corrupted aggregate fixtures were rejected after updating
their result/CSV/run-log hashes consistently, demonstrating numeric/scope guards
rather than trivial checksum rejection. Only the mutation fixtures skip chart
checks because they deliberately change the input result; the default CLI and
both real-output audits require the chart files/provenance. No original output,
source or repository file was changed. Logs labelled Python312/313 each show27
tests and OK; plain saved logs alone do not authenticate the interpreter binary
or remote CI. The reviewer also previously executed27manufactured unit tests.

Actual output SHA256:

- RESULTS_SUMMARY.json: `dd5e688effd6324fbc07315722c3c3603ca44ae25a6f535ece79b4e0a7ce5846`
- REFERENCE_TABLE.csv: `5234470d6af5c0e3c449f8d8de0762c6e496494c4c5b5b02d8be45650e960f66`

No critical finding remains for the provisional local-derived period-reference
scope. Original UN RDA scale/bytes, a current-calendar or personal forecast,
disease-free intensity, biomarker adjustment and updated SHFM5year prediction
remain explicitly **unverified/unperformed**, not established by this audit.
