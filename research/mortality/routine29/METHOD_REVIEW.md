# Independent pre-result method/code review29

2026-10-29. Read-only review of protocol7800458, LOCK, renal/BMI metadata,
eGFR rationale, features.py, cohort.py, retrieve.py and routine29/model.py.
No cohort records, comparative results, model execution, data retrieval, ZIP,
account, agreement or contact used. Numerical module was authored by this
reviewer, so its independent validation rests on explicit synthetic tests;
this is not falsely presented as a separate third-party audit of that module.

## Scope and current verdict

The prespecified question is three proteins beyond a compact renal/routine
core in the REUSED NHANESIII cohort. It is not an independent validation.
No scientific inconsistency found in the completed pieces reviewed so far.
At this initial review checkpoint reader.py and analyze.py were not present;
these remain pending implementation inspection before the final verdict.

## Prespecified contrasts and measurements

Primary routine_all3 minus routine, same complete-case domain, shared survey
weights and five whole-PSU outer labels. Simple/simple_all3 are re-estimated
on this SAME routine-common domain; saved28 scores on3857 are not the proper
cross-domain comparator. No individual-protein model selection or pruning.

Eleven routine terms correctly represent six previous clinical terms plus
BMI,HbA1c,eGFRlow/eGFRhigh,log2UACR. Calibrated serum creatinine is not also
included as a second predictor. Diabetes diagnosis is not duplicated with
HbA1c; cystatin is not used in the eGFR comparator.

CDC serum calibration .960*CEP−.184 acts on serum mg/dL only. Nonpositive
corrected results remain missing.2021 adultcreatinine-only eGFR formula in
features.py matches the metadata rationale (142,min/max,sex-specific k/alpha,
−1.2, .9938^age, female1.012), with no race multiplier. It is an estimate,
not measured GFR or proof of accurate historical-individual renal function.
A fixed60hinge is recorded before results, not selected by comparisons.

UACR100*UBP/URP uses urine units ug/mL and mg/dL. Primary requiresURP>=10,
which excludes the documented statistically-suspect7.9 replacement before
outcomes; serum calibration is not applied to urine. ReleasedUBP.4 is retained
and flagged, not declared aboveLOQ or imputed anew. PrimaryBMI retains valid
reported/edited values; anthropometric substitutions and GHP method are
reported. Unknown/missing source codes stay missing. Strictmetadata mask
BMIflags0/Hbmethod1–3/UBP>.4 is secondary only and does not certifyQC.

## Source and cohort boundaries

Study hashes and original public qualification gates run before participant
parsing. Newmanifest is bound to the lock and old source manifest; all files
are checked for byte lengths/SHA, names and official release-index hrefs.
Original1a LAB/EXAM is distinguished from3a repeats. The larger EXAM cap220MB
is justified by documented31311*6235 nominal bytes, not an unbounded download.

Oldclinical/protein eligible3899 and oldcommon3857 are preserved and checked.
New routine mask includes the six oldclinical terms, so all primary arms
share a full common domain without silently normalizing missing predictors.
Source join is one-to-one and checks survey-ID universe. Different marker
coverage/population must be disclosed; complete-case selection remains a
limitation, not fixed by the source's original survey weights.

## Numerical learning and uncertainty

New model.py reuses immutable binary28 numerical primitives without mutating
its registry or files. Simple arms use exactly first6features; routine arms
use all11; rank/tail transforms are trained only inside each actual split.
The new routine module currently has14synthetic tests, including
holdout outcome/feature leakage checks, PSUs, zero weights, constant columns,
source input mutation and preserved binary28 registries. Sameinner partition
simple-model predictions match the preserved algorithm to floatingprecision.

Nested training chooses lambda1/10/100 by aggregate weighted innerOOFBrier,
ties within1e−12 favor stronger ridge. Inner salt routine29:inner:{fold};
normalization retains raw training n including zero-weight rows as specified.
All preprocessing/inner tuning/outer fitting must be repeated inside500paired
rescaled-PSU draws. Fixed partitions remain conditional design choices.

Five additional outcome-blind partition studies are correctly prespecified
as dependent sensitivity ranges/directions, not pooledbootstrap intervals
or independent validations. Strictmetadata analysis must fit/evaluate both
primary arms on its own paired subset with the original outerlabels, not
compare its raw scores to the main population.

## Limits that must remain in the report

Known positive28 result, reuse of same cohort/programme history; reported
protein QC/LOQ/calibration/storage limits; source substitutions and urinary
nondetect handling; complete-case selection; retrospective eGFR assumption;
finite500bootstrap MonteCarloprecision; conditional fixedpartition inference;
no clinical usefulness, external transport or universallyoptimalpanel claim.

## Completed reader/analyze inspection (pre-results)

The now-present reader.py streams source text, keeps literal decimal values,
checks exact INPUT spans independently of LENGTH/FORMAT, validates nominal
source width/MISSOVER and official row counts, then requires exact equality
with pandas.read_fwf. Selected end1970 for LAB and1533 for EXAM includes the
requested metadata flags. All actual raw source bytes/newlines are hashed;
unknown tokens, missing required ranges and duplicate/invalid IDs stop import.
Nothing executes SAS or silently divides numeric data by a FORMAT decimal.

analyze.py confirms all six point-estimate arms use the same core/protein/y/w
arrays and folds from eligible3899 before the common mask. Primary routine
pair is refitted after all preprocessing/tuning in each of500 shared-weight
replicates. Selection counters are updated only after successful paired fit
and scoring; per-metric undefined diagnostics remain separately reportable.

Each of five extra outer allocations is a new outcome-blind hash over eligible
clusters, followed by the same common restriction and nested procedure for
both routine arms. All records/failures remain in partition_sensitivity;
ranges are separate and partition_sensitivity_pooled_into_CI=False.
These are dependent sensitivity studies, not a new nominal-confidence band.

Strict metadata conditions are intersected with the primary common domain.
Both routine arms are refit and evaluated on exactly that subset, with original
outer labels. Missing fold/event support causes explicit undefined status,
not a favorable reallocation. Its paired contrast is not a cross-population
comparison of absolute Brier/AUC. No bootstrap or inferential winner is added.

Point calibration bins/parameters and all500primary delta rows are exported
through marked aggregate stdout, enabling independent arithmetic/percentile
checks without ZIP or another model calculation. No participant predictions,
IDs or empirical ECDF support values are exported by the model records.

### Nonblocking reporting recommendation sent to parent

The plan requests exclusion-reason overlap. Current missing_routine fields
count11predictor columns, including duplicate eGFRlow/high missingness, and
combined removed-from28 count, but not an exclusion-pattern table. Add grouped
clinical/protein/BMI/HbA1c/eGFR/UACR overlap counts (and common-domain flags)
so the new population can be understood without summing overlapping columns.
This requires only aggregate reporting, not a new predictor or scientific
selection change. The parent was notified before comparative results.

## Final pre-result verdict

No blocking method/implementation error found. Reader, common-case comparison,
whole-PSU inner/outer learning, paired full refitting, partition sensitivity
and within-subdomain strict metadata comparison agree with the locked plan.
The study remains exploratory on reused data and must retain assay/eGFR/
missingness/partition limitations. No source records or comparative results
were read by this review; no files in the repository were changed here.

## Parent resolution before comparative outputs

The requested aggregate reporting was added before reading model comparisons:
exclusion_reason_counts, exclusion_reason_patterns, patterns within old28common,
common-domain source flags and UBP0.4 counts. eGFR low/high are grouped as one
measurement-level exclusion reason. No predictors, selection rules or primary
contrast changed. The numerical module's author explicitly disclosed the scope
of self-review; independent result/CSV arithmetic is a separate audit.
