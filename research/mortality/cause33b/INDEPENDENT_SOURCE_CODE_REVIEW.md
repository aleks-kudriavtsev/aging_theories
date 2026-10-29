# Cause33B: independent source-contract and shared-draw code review

Reviewed2026-10-29 at **science `ab09ca0b605b44544f49d21d1161c61609ef92d1`**, protocol **`d772b874198a9f7739225bbce5d1b6b6dfcd0c08`**. Scope: actual `cause33b/cohort.py`, `analyze.py`, `retrieve.py`, CONTRACT/LOCK and inherited source/model/bootstrap semantics. No33Bresult contents, live logs, CI polling, rawrecords, sourceGET/ZIP or refits. Only16pure projection/comparator checks on synthetic dictionaries were executed. No repository edits; `/tmp/verify33B.py`, laboratory report and presentation untouched.

**Decision:** no remaining scientific blocker found in the reviewed source projection, cause construction or shared-draw orchestration. This is pre-result code evidence, not a claim that actual runtime source qualification, model fits, intervals or final-headCI passed.

## Protocol, source and code pins

Independently hashed all **20metadata pins** locally; allmatch. All **99historical scientific Python hashes** remain unchanged;4new cause33b Python modules give **103total**. Source intake remains inherited31:50required inputs,25historical data fingerprints and two separately recorded prior-denied/not-retried PDFs. New33Bwrapper binds its own schema, protocolcommit and LOCKSHA to exactbase-manifestSHA and filefingerprints. `verify_wrapper` runs before the old parser and again after source/endpoint qualification; loadedbasemanifest must equal verifiedbase. No new source or measurement policy is introduced.

Actual33Bsource contract reruns the preserved33Agate on the loaded eligible/common frame, checks source qualification, complete early causes, permitted comparisons, absence of blockers and supported targets in exactorder **`['002','001']`**, then compares all retained scientific aggregates against the frozen33Aresult. The frozenfile itself is LOCK-pinned. This verifies intended domain/count/fold reproduction before fitting, rather than accepting a fresh self-consistent manifest as a new source identity.

The runtime source check still depends on the preserved31raw-byte/schema/marginal/reader guards. I checked local code and metadata bytes only, not current source files or any new computed gate result.

## Exact scientific projection versus runner envelope

Livegate must have **exactly19scientific top-level fields**, schema`cause33-count-gate/1`. Frozen33Aoutput must have those19plus exactly9declared runnerfields, schema`cause33-qualified-counts/1`. Unknown or missing top-level fields fail; the frozen schema alone is normalized for comparison. Kept fields:

`schema, gate, endpoint, eligible, common, per_cycle, causes, unknown_cause, defined_cause_ranking, priority_top5, proceed_codes, accounting, design, folds, weight_rule, fitted_models, predictive_metrics_computed, participant_records_exported, threshold_scope`.

Removed frozenrunner fields are only `plan_commit, lock_sha256, runtime_python, prior_cohort_use, comparison_models_fitted, independent_binary_probabilities_summed, joint_integration_and_Fine_Gray, historical_scientific_files_verified, scientific_python_sha256`. Numeric cause counts, priorities, support failures and fold fingerprints are **not** stripped. Source-fingerprint identity is additionally bound by the wrapper/inheritedloader, separate from this projection.

The recursive comparator requires **identical PythonJSON types first**. Dictkeysets, listlength/order, integers, booleans, strings, fingerprints andnull are exact. Therefore1≠1.0 bytype andTrue≠1; causecodes/counts cannot exploit numeric tolerance. Nonfinite or unsupported non-JSON values fail. Only float-versus-float derived reductions use `math.isclose` with locked `rtol=1e−12`, `atol=1e−9`:

`abs(a−b) <= max(1e−9, 1e−12 * max(abs(a),abs(b)))`.

This is **not** numpy's additive/asymmetric tolerance formula, and not a rule that every absolute difference must be≤1e−9. Largeweightedtotals can legitimately have a larger allowedabsolute rounding difference. Counts/types/sourcehashes remain strict. The CONTRACT's tolerance fields must themselves be floattyped and equal the fixed constants, and the CONTRACThash is pinned.

Proof exports floatchecked/different counts, maximumabsdifference, comparedfieldlist and **both** canonical scientific hashes (`sort_keys=True`, compactseparators, allow_nan=False). These hashes may differ when tolerated roundoff is accepted; do not require equality as an extra unregistered gate. Identical-valued signedzero can also serialize differently. The comparator, not hashidentity alone, establishes the declared scientific equality.

Syntheticchecks **16/16passed** for schema/projection, knownrunner stripping, unknown/dropped fields, exactcount/boolean types, listorder, label/field changes, nonfinitefloats, and independent mathematical tolerance-boundary cases atsmall/large scales. This tested AST-extracted pure functions, not source loading or model calculations. There are30synthetic testdefinitions in repository cause33b tests; I did not run that suite and make no CIclaim.

## Actual target labels and controls

`cause_targets` constructs **new cause-specific labels**, avoiding inherited31's allcausey helper. It requires numeric nonboolean status/cause/timecolumns, exactboolean commonmask/shape/index, status0/1, finite nonnegativetime, integralknowncauses1..10, no cause foralive and no alivefollowup<96. Earlycommonunknowncause blocks. Returned arrays are `dead==1 & t_exam<=96 & cause==k`, in unchanged common-roworder, for002then001. Unknowncause after96does not create an earlyunknownlabel.

The locked prior countchecks are N6044, eligible6546, cancer002events176/controls5868, heart001events148/controls5896; the33Agate separately reproduces allcause586. These are committed prior aggregate expectations, not newly read stageBresults here. Controls include known competing earlydeaths, laterdeaths and observed96monthsurvivors. No cause-specific cohort restriction, survivor-only AUC, censoring of competingdeath, or unknown-to-zero substitution was found.

Run checks exactN/event/controlcounts before fitting. Numericalvalidation requires binarytargets, samealignedX/proteins/weights/rows/clusters/folds and mutuallyexclusive indicators. It verifies commonweights equal the eligible-frame weights selected by thatmask, clusteralignment, equaleligible/commonPSUsets and exact31outermapping. Source-measurement policies, W6, proteinorder and transformations remain inherited/unmodified.

## One shared rescaled PSU draw

`analyze_arrays` creates one RNG with lockedseed20261105. Every attemptedreplica calls `design_bootstrap_weights(eligible,rng)` **exactlyonce before either cause loop**, selects the onecommon `wr`, computes its sharedsummary/fingerprints, then sends thatsamewr to everycause and everyarm. Zero-weight rows and originalfoldIDs remain. There is no conditionalredraw/reseed or separately sampled cause/armdraw.

Inherited sampler uses sortedstrata and eachstratum's actual publicPSUinventory. Drawm−1unitswithreplacement, multiplier`count*m/(m−1)`. For everystratum `sum(count)=m−1` and summedPSUmultipliers=m. TwoPSUstrata give0/2; threePSUstrata allow0/1.5/3. Weightedrowtotals neednot equal pointweighttotals. NativeW6 evaluationweights remain the shareddrawweights; training normalization is inherited fit-local and allocates a new weightarray, rather than mutating sharedwr.

Exports include SHA256 of exact alignedweightbytes `np.asarray(w,dtype='<f8').tobytes(order='C')` for eligibleandcommon vectors, plus sharedsumw/sumw²/positiveweightcounts/PSUs and per-cause weighted/positiveevent/noncasecounts. **No weightvector or participantidentities is exported.** A single fingerprint perrep plus reviewed one-array orchestration establishes the available evidence; without rawvectors/publicPSUdrawcounts, the aggregateauditor cannot independently regenerate that digest. Do not describe it as independent rawweight reanalysis.

Samearm/outerpreprocessing metadata are compared exactly across allavailable matchingcausefits atpoint and everyreplica. The result explicitly limits this check to availablepairs after fitfailures; missingpairs are not claimed verified. The check covers means/scales/constants/featureindices/names and geometry/weightaggregates, not coefficients or selectedλ, which may correctly differ byendpoint.

## Failure/penalty-count nuance for the independent auditor

Cause/arm fitting errors are local; metricerrors are separately prefixed `fit|` or `metric|`. ValidBrier/undefinedAUC or calibration remains represented by metric-localnulls/statuses; pairedCSV has individualBrier/AUCstatus. All500attempts remain; marginal exploratory bounds use finitevalues permetric without FWE/interactionclaims. No probabilities are summed.

**Parameter rows are appended immediately after successfulnestedtraining, before scoringmetrics.** Therefore bootstrap penaltyselectioncounts track completedtraining, even if subsequentmetric evaluation fails. Under this actualcode:

`sum(selected_penalty_counts[cause][arm]) = 5 * (500 − number_of_fit-prefixed_arm_failures)`.

It neednot equal5×computedmodelstatuscount when metric-prefixedfailures exist. Pointparameters can similarly exist for a metric-undefinedarm. `_check_shared_preprocessing` can inspect those successfultrainingrows. The separatefailure summaries and hash/CSVoutputs make this distinction auditable.

## Reviewed final-science hashes

- `cause33b/cohort.py`: `a723340edaf441cb2eeb02a0d6ffe24bbeb30a548c809bf3fce7a1f0c65af673`
- `cause33b/analyze.py`: `79d82765db720afd774a45c93750fe3e939282f88487154addd899e104d96f17`
- `cause33b/retrieve.py`: `2b08cf53adfa7eb19052fab8c31d4e11f538027dd19aa66badef06e0978032b8`

Remaining verification belongs to the parent's final-headaggregate/CIaudit: actual50source/runtimefingerprints, projectedfloatdifferences, models, metric-localfailures,500pairedpercentiles and exactstdout/filehashes. This review does not change the preregistered tolerance, policies, models or deliverables.
