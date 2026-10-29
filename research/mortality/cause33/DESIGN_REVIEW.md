# Cause33 design review: separate8-year recorded-cause binary comparisons

2026-10-29. READ-ONLY design for the latest user priority. No repository edit,
retrieval, participant-file read, outcome count computation, fit or model metric
was performed. Joint competing-hazard integration and Fine–Gray are DEFERRED;
this note does not authorize either. It supersedes the earlier direction toward
a joint/Fine–Gray benchmark. Existing33-search notes are historical notes only.

## Recommended sequence: two distinct commits and execution stages

**StageA now: source/count/support gate ONLY.** Commit its schema, cause-code
policy, candidate ranking, thresholds, fixed folds, tie rules and stop behavior
BEFORE reading new cause counts. Reuse source-qualified SSCARD1999–2004 plus
public2019 LMF and the EXACT31 primary common domain. Produce observed source
marginals, all candidate counts/burden/support, deterministic ranking and up to5
qualified causes. Do not compute AUC, Brier, calibration, model coefficients,
penalty grids or predictions in this stage. A support table is not model evidence.

**StageB later: comparative plan AFTER StageA.** Save the StageA commit/output
hash and explicit selected codes/reasons, then commit the model/contrast/uncertainty
plan BEFORE fitting anything. State that cause selection uses full-cohort outcome
counts/burden and the reused cohort; neither selection nor subsequent CV is an
untouched/independent validation. Do not call the later plan outcome-blind.
If no cause passes, StageA ends with a concrete feasibility result, not invented
models, changed thresholds, a favourable rehash or forced inclusion of five causes.

## Source, endpoint, assay and design gate

Reuse exact31 source fingerprints and code/metadata versions, not a new guessed
SSCARD URL or new weights. Relevant existing code:
- b2m31/cohort.py: verify_sources, validate_sscard, six_year_weights, domain_masks,
  validate_design_inventory and load; it preserves the underlying-cause column
  because the full public mortality table is merged, although31 used all-cause y.
- validate_temporal08.py: independent_mortality_reader,1-based columns eligibility15,
  vital status16, UCOD_LEADING17–19, interview43–45 and exam followup46–48.
  Its second-reader audit is parsing agreement, not proof of clinical definitions.
- b2m31/features.py: unchanged11-term routine core and verified cycle-specific
  creatinine mapping; proteins ordered SSCYST(Cys),SSB2M, log2 reported mg/L.
- b2m31/model.py and analyze.py: existing pure binary/nested/PSU primitives.

Check these items before returning any supported cause:
1. All3 cycles are1999/2001/2003 (baseline1999–2004), eligible public2019 linkage,
   known binary vital status, screening age40–75, MEC examination, valid sex,
   documented masked design and positive CORRECTED stored-serum weights.
2. Verify that these source years expose the full UCOD_LEADING001..010 partition.
   Never infer availability from newer2015+ files, which have coarser public causes.
   Record the codebook version/URL and current retrieval status;31's two denied
   mortality PDFs remain explicitly unavailable under its amendment, not newly
   verified. Do not retry denied documentation without relevant new authorization.
   Existing reviewed code/source mappings can be cited honestly; an unresolved
   required code meaning is a gate blocker, not an invented synonym.
3. Nonmissing underlying causes on deaths must be supported integer codes1..10;
   reject malformed/unsupported values. Survivors must not have a death cause.
   Underlying cause is ONE mutually exclusive UCOD code. Do not substitute the
   separate diabetes/hypertension contributing-cause flags as underlying endpoints.
4. Time origin is exam, tau=96months inclusive. Alive before96months means STOP
   for this uncensored binary design; do not exclude those people, shorten tau or
   silently add IPCW. Death after96months is a known target noncase at96months;
   death by96months from another cause is also0, not censored. Month0 deaths would
   be ordinary binary events, without artificial time shifting.
5. Preserve the exact31 common domain for EVERY cause and arm: finite11 routine
   terms, both positive finite proteins and both detection flags0. No cause-driven
   missingness changes, LLOD sensitivity, subgroup restriction or assay substitution.
   Flags0 indicate the released measurable range, not individual QC/modern assay
   traceability. Keep storage/assay timing and prior cohort-use limitations.
6. Use pooled six-year weights as31, NOT32's native-period split weights:
   w=(2/3)*WTSSCB4Y for1999/2001, w=(1/3)*WTSSCB2Y for2003. Late W4 missingness is
   expected and must not remove2003 participants. Preserve raw W2/W4/W6 provenance.
7. Cluster=(SDMVSTRA,SDMVPSU), not PSU alone;87 masked PSUs across43 strata were
   recorded by31, including one3-PSU stratum. Require the eligible survey design
   inventory and common-domain inventory to agree, or explicitly block the planned
   eligible-derived inner-fold equivalence. No singleton or dropped-PSU repair.
   Future resampling draws the full qualified eligible frame BEFORE common masking.

Known31 aggregates may inform feasibility planning only: n6044,586all-cause deaths
by8years,87PSUs/43strata, common-cycle n1889/2131/2024. They do NOT establish cause-
specific counts or support. Do not infer those counts from marker associations.

## Endpoint policy: what each binary target means

Candidate codes are001..009, after primary source definition verification:
heart diseases, malignant neoplasms, chronic lower respiratory diseases,
unintentional injuries, cerebrovascular diseases, Alzheimer disease, diabetes,
influenza/pneumonia, nephritis/nephrotic syndrome/nephrosis. Use exact release
labels in the final registry. These are coarse recorded groups: heart is not IHD,
respiratory is not COPD alone, Alzheimer is not all dementias and renal is not
CKD alone.004 and008 must not be silently excluded if the user's specified
candidate range is001..009; label injuries/infection plainly rather than calling
all nine noninfectious diseases.

For target k: y_k=1(dead=1, exam_followup<=96, UCOD_LEADING=k); all other known
outcomes are0. Code010 is heterogeneous and excluded from CLINICAL RANKING, but
its deaths stay as competing noncases in every target model and all-cause totals.
Do not drop010 participants or assign010 to a marker/disease module.

Count missing underlying cause among known deaths separately from coded010.
Never drop/censor such deaths or guess their disease. If they exist, either:
(a) explicitly target RECORDED-CODE cause risk, with uncoded deaths retained as
noncases and the missing-cause limitation stated; or (b) require complete death
cause ascertainment for a stronger underlying-cause estimand and block that
claim if unmet. Choose/pin the policy before counts; do not silently imply that
coding an unknown death0 proves the person's true cause was not k.

Each fitted binary probability estimates P(recorded death from k by8years|X)
in this received complete/measurable domain. It is a fixed-horizon CIF comparison,
not a cause-specific instantaneous hazard model. Competing deaths are observed
noncases. Separate models must NEVER be summed into all-cause risk or survival,
or presented as a coherent joint state vector. Shared markers may enter every
cause model; there is no exclusive B2M→renal/Cys→cardiac assignment or causal claim.

## StageA ranking and conservative computational support

Show ALL9 candidate rows, even excluded ones, plus010/unknown as unranked source
accounting rows. Let E_ik be the target indicator above. On the SAME common domain:
- raw events=sum(E_ik);
- weighted8-year burden=sum(w_i E_ik)/sum(w_i);
- share of weighted8-year mortality=sum(w_i E_ik)/sum(w_i D_i8);
- event-bearing PSU and stratum counts, nonevent support, and event-weight Kish n
  ((sum wE)^2/sum(wE)^2), largest event-weight contribution as diagnostics.
These are received-domain estimates, not unqualified national disease rankings.

Rank by descending weighted8-year burden; ties by descending raw events, then
ascending code. Select the first at most5 candidates PASSING the fixed support
gate. Publish the full raw ranking and skipped reasons. Never rank by AUC,
Brier improvement, p-value or a previously favourable marker. Report fewer than5
when support permits fewer; failure is an informative result.

Recommended conservative thresholds to pin before counts (root may set a
specified alternative before execution):
- total >=50 positive-weight target-event rows in >=15 event-bearing PSUs;
- at least20 positive-weight nonevents overall and both classes in every fit;
- each fixed outer scoring fold >=2 target events in >=2 event-bearing PSUs;
- each outer training set >=30 target events in >=10 event-bearing PSUs;
- every actual inner scoring fold >=1 target event in >=1 event-bearing PSU;
- every actual inner training set >=20 target events in >=6 event-bearing PSUs,
  and >=20 positive-weight nonevents. Check all15 actual inner train/score pairs.
Also expose event strata and weight concentration even when the thresholds pass.

These are conservative COMPUTATIONAL/REPORTING thresholds, not a power analysis,
EPV approval, calibration adequacy guarantee or clinical readiness. Penalized13-
term models can remain very uncertain with50events. PSU bootstrap can still
produce undefined fits/AUC; no threshold makes every replicate valid. If root
prefers a lower floor, label its evidence appropriately rather than implying
that>=10events validates a clinical model.

Assign a SINGLE outcome-independent whole-PSU5outer×3inner partition by pinned
hash salts BEFORE common restriction/count-based support checks. Recommended new
salts cause33:outer:v1 and cause33:inner:{outer_fold}:v1, or explicitly reuse31's
existing salts; choose once before counts, never compare/reselect allocations.
Check eligible/common mapping equivalence. Do not move a PSU or change folds
because a candidate narrowly fails the support gate.

## StageB compact comparative design, not yet authorized to execute

For every selected cause independently, use unchanged11-core and log2(Cys,B2M)
with all four31 arms: routine, routine_b2m, routine_cys, routine_cys_b2m. Same
participants/folds/weights in all four and all selected causes. Primary per-cause
contrast: routine_b2m minus routine, weighted8-year Brier. Prespecified secondary:
routine_cys_b2m minus routine_cys. Report all four arm metrics and both contrasts;
no best-protein selection or exclusive disease assignment. AUC and unchanged-risk
calibration are secondary diagnostics; rare-event undefined counts are mandatory.

Five outer whole-PSU folds,3inner per outer; training-only weighted mean/SD,
constant scale1, no imputation/clipping/new features. Lambda1/10/100 by pooled
weighted inner-heldout Brier, ties1e-12→larger lambda. Fresh outer training refit.
Use a NEW cause33 wrapper/artifacts without changing31/32 or frozen models.
The existing31 model accepts an arbitrary binary y, so its tested numerical
primitives are suitable; avoid mutable global-registry/salt monkeypatches.

500 paired rescaled eligible-frame PSU draws shared by all selected causes/arms;
relearn all preprocessing/tuning/outer fits. Keep original rows/PSU/fold labels
with zero replicate weights. Retain every attempt, arm-local/metric-local errors
and primary/secondary validity separately. No redraw to obtain500 successes.
Intervals are conditional on selected endpoints, fixed partition/forms, reported
measurements and prior analytic history; no external/untouched claim. Per-cause
95% intervals are not simultaneous coverage across up to5causes/two contrasts;
either keep all results explicitly exploratory or predeclare multiplicity
control in StageB. Calibration diagnostics never update predictions.

## Practical tests and aggregate artifacts

StageA tests, before any model code is invoked:
- every allowed/invalid/blank cause code; dead/alive inconsistency; contributing
  diabetes/hypertension flags cannot change underlying-target labels;
- horizon boundary95/96/97months, competing death0,010/uncoded policy, survivors
  <96blocked, no dropped deaths/no cause-based denominator change;
- correct historical W6 formula with missing lateW4; native32 weights rejected;
- fixed common mask across all9codes; unchanged identifiers and source checks;
- full strata/PSU inventory,3-PSU stratum, whole-unit nested folds, exact support
  counts, failure of a single inner-training/scoring cell, no rehash after failure;
- weighted-burden ranking opposed to raw-event ranking, deterministic ties,
  unsupported high-burden cause skipped with reason,010 never selected,0..5causes;
- eventPSU threshold distinguishes clustered events from identical event totals;
- gate stage imports/calls no fitting, metric or prediction routines (mock them
  to raise); deterministic JSON/CSV hashes and aggregate-only outputs.

StageA outputs: source manifest/codebook availability, observed release marginals,
cohort/endpoint flow, all-code burden/support table, fixed-fold support aggregates,
selected codes/ranks/reasons and exact source/code/protocol hashes. No rows/IDs/
prediction vectors. StageB adds separately committed analysis lock, aggregate
model parameters/OOF metrics/calibration bins and paired bootstrap CSV, plus
independent arithmetic/hash audit and scientific CI. No accounts, agreements,
restricted data, holder contacts, merges or deployments are implied.

## Final pre-count addendum —2026-10-29

Reviewed the now-written cause33/GATE_PLAN.md and SUPPORT_RULES.json before
protocol commit and before new cause counts. **No scientific design blocker
found in these final choices.** They take precedence over the alternatives in
my earlier recommendations above:

- The final gate ranks all defined001..009 first, takes the top5 PRIORITY
  candidates, and then reports their support failures. It does not backfill with
  rank6 or below. This differs from my proposed first5 supported candidates but
  is a legitimate, clearly stated question about the leading burden groups.
  Report fewer than5 without implying that a failed candidate has zero risk or
  that a lower-ranked cause could not support another separately scoped study.
- The final conservative JSON minima are50 event/50 noncase rows overall,
  >=20 event-bearing masked PSUs, outer evaluation>=5 of each class, outer
  training>=20 of each class, inner evaluation>=3 of each class and inner
  training>=10 of each class. These differ from my suggested count/PSU limits;
  none of these alternative pragmatic minima is a universal optimum or a power,
  EPV, calibration or clinical-validity certificate. I recommend no threshold
  change merely to conform to my earlier arbitrary values. Retaining precise
  failed cells and later bootstrap undefined counts remains essential.
- The final design reuses31's fixed eligible-derived b2m31:outer and
  b2m31:inner:{outer} salts rather than inventing new salts. This preserves the
  comparison history and is acceptable when explicitly described as internal
  work on reused data. The eligible/common inventory equivalence gate prevents
  silently changing the inner maps; no rehash or fold repair after counts.
- The final plan chooses the stronger underlying-cause ascertainment policy:
  ANY missing UCOD among common-domain deaths at/before96months blocks later
  cause comparison, while remaining visible in accounting. This resolves the
  earlier note's two possible policies without silently treating unknown early
  causes as negatives. Missing cause after96months does not change the8-year
  target. Coded010 remains an observed competing noncase but is not ranked as a
  clinically directed cause.

StageA and StageB remain correctly separated. The earlier phrase “StageB not yet
authorized to execute” refers ONLY to the methodological requirement to freeze
a separate pre-comparison plan after source/count qualification; it is **not a
request for new user permission**. The user already authorized research work.
Proceed within that authorization once the scientific gate and committed plan
are satisfied. Describe the Git plan accurately; no claim of external protocol
registration or unchanged outcome-blind preregistration is implied.

Suggested diagnostics need not become additional gates: per-cell event-bearing
PSU/stratum counts, event-weight Kish n and maximum event-weight contribution
help readers see concentration even when row thresholds pass. Source-only gate
tests should assert that fitting, tuning, metric and prediction APIs are never
called; top5-before-support/no-backfill and the final JSON thresholds need exact
boundary tests. Weighted event total versus weighted incidence gives identical
rank order here because every candidate uses the same fixed denominator.

No cause counts, ranking values, comparative metrics, fits or retrieval were
performed during this review. No repository file or earlier model was changed.
Joint integration and Fine–Gray remain deferred.
