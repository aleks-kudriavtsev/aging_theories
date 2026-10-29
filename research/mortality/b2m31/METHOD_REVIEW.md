# Independent pre-result scientific/method review31

2026-10-29. Read-only review of b2m31/ANALYSIS_PLAN.md, features.py, model.py,
analyze.py and /tmp/b2m31_source_weights.md. No participant files, existing/new
comparison results or metrics were opened; no numerical analysis executed.
The protocol is now committed6774fbe. The cohort loader was not yet present at
this review: its source, flags, full-design and endpoint checks remain required
before the first empirical comparison rather than being certified here.

## Verdict

The proposed reused-cohort exploratory study is scientifically interpretable
within its declared scope: a new8-year weighted binary ridge model estimates
whether reported B2M adds information to the fixed11-term routine core, with
an explicitly secondary assessment beyond cystatin C. It is neither a new
untouched cohort nor external validation of freeze30 or of the original engine.
No critical leakage or statistical-model inconsistency found in the reviewed
feature/model/analysis code, subject to the source/design gates below.

One concrete pre-result inconsistency was identified and corrected: the initial
analyze.py seed20261103 differed from plan20261102. Read-back now confirms
SEED20261102 and PLAN_COMMIT6774fbe. The protocol and interpretation now also
exclude flag2 from the reported-values sensitivity, as discussed before results.

## Source/assay qualifications independently checked against official metadata

Current [SSCARD_A](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.htm)
is the July2022 revised release. It removed erroneously included records and
recreated its specific weights without changing retained measurements. Both
proteins are serum mg/L measured2018–2020 by Siemens Dimension Vista1500, from
biological baseline1999–2004 specimens. Flag0 means measurable;1 below LLOD
uses LLOD/sqrt2. Upper-LOD flag2 replacement semantics are not described, so
excluding2 from the limited sensitivity avoids pretending an unknown censored
representation is an exact concentration. Detection flags are not LOQ or
individual-QC certification; no older cystatin calibration should be imported.

The same codebook's non-pristine warning concerns glycated-albumin performance;
it is not evidence that all non-pristine Cys/B2M are invalid. Keeping both
pristine categories in the primary domain while reporting SSPRIS is defensible.
No unplanned pristine-only subgroup or outcome-selected cutoff is needed.

Current [LAB18 analytic notes](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/LAB18.htm)
explicitly give StandardScr=1.013*original+0.147 for1999–2000. features.py uses
this only for1999 and does not reuse NHANESIII's different equation or correct
urinary creatinine. The2001–2004 column/scale claims still require the loader's
explicit source-field/unit verification; names must not be guessed.

## Pooled weights and PSU design

The official [NCHS weighting tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx)
combines1999–2004 with2/3 of the supplied1999–2002 four-year weight plus1/3
of2003–2004 two-year weight. Its footnote applies the same structure to subsample
weights. Applying that documented rule to revised WTSSCB4Y/WTSSCB2Y is justified;
WTSSCB2Y/3 for all cycles would be wrong. Separate within-cycle normalization
before pooled scoring would also change the target and must not occur.
The model's training-only likelihood normalization is a distinct prespecified
regularization convention and does not replace raw survey weights in metrics.

Clusters must be unique(stratum,PSU) pairs. SDMVPSU1/2/3 is permissible in the
first cycle; no hardcoded2-PSU assumption belongs in31. The new bootstrap routine
uses each stratum's actual m clusters and draws m−1 with factor m/(m−1), rather
than borrowing a two-PSU-only implementation. Whole copies retain fold identity.

Crucial qualification: compare the eligible precomplete pair inventory with
ALL matched-cycle DEMO pairs BEFORE any primary comparison. Stop if any pair
is lost or a stratum is singleton; no relabelled strata or arbitrary singleton
rescue. The final committed protocol includes this gate. The source loader must
implement and export actual full/eligible inventories and equality results.
This removes missing-PSU ambiguity; it still does not make this nonlinear,
training-and-tuning bootstrap an official NCHS variance estimator. Describe it
as a conditional rescaled masked-PSU refit bootstrap, retaining assay/missingness,
prior-selection and partition limitations. Survey weights do not correct extra
core/protein/flag/linked-eligibility complete-case selection.

## Domain, features and nested evaluation

- All FOUR primary arms require the same11 finite routine features and both
  proteins positive finite with both flags0. Requiring cystatin for the B2M
  contrast is a deliberate common-domain restriction, not an estimate for all
  persons with B2M. Report weighted coverage and cycle/domain changes.
- Smoking, BP mean, BMI/HbA1c and renal transformations are source-specific and
  separate from NHANESIII. Missing/refused smoking codes are not normal/never;
  BPXSY4 is not silently substituted. There is no added raw-creatinine term or
  cystatin-derived eGFR in the comparator. Thus the incremental-B2M question
  controls the specified creatinine-eGFR/UACR core without claiming a complete
  universally adequate clinical panel.
- Screening age and original exam followup are not serum assay/release dates.
  All eligible survivors must support96months; missing endpoint or early censor
  cannot become a nonevent. Source-specific routine released replacements remain
  labelled; no NHANESIII urine-creatinine floor or HbA1c crossover is imported.
- Both outer and inner folds are whole-PSU, deterministic and shared across arms;
  feature means/scales, lambda search and outer fits use training rows only.
  The registry selects only intended protein columns and never mutates earlier
  model registries. No old prediction or coefficient enters these fits.
- Lambda selection is pooled weighted inner-heldout Brier with the prespecified
  ties; fixed-partition outer predictions are evaluated once on each heldout row.
  Retained zero-weight rows in resamples keep partition identities stable, while
  weighted moments exclude their influence. Training normalization including
  zero-weight rows matches the explicit objective convention.

## Uncertainty, contrasts and detection sensitivity

500 paired bootstrap draws use the same qualified eligible frame and draw
weights in allfour arms; each reruns preprocessing, tuning and fitting. Primary
B2M-minus-routine and secondary B2M-beyond-Cys deltas use paired predictions and
weights. Per-model, per-contrast and metric-local failures are retained: an
undefined secondary arm need not erase a defined primary Brier contrast. AUC
or calibration failures are not silently treated as successful estimates.
Calibration diagnostics never change scored probabilities. Fixed-fold refitting
intervals are conditional on the selected hypothesis, model forms, partition,
reported measurements and programme history; they do not remove post-selection
or calibration/storage uncertainty.500tail quantiles have Monte Carlo limits.

Primary and secondary must remain distinct. Both intervals may be reported,
but neither an additional secondary success nor the collection of four models
establishes multiplicity-free confirmatory significance or a universally best
panel. B2M was selected after earlier programme results; that is already
properly acknowledged in the protocol.

The sensitivity now uses complete core and known flags0/1 on BOTH proteins,
excluding2, with shared participants across its four arms and point estimates
only. It is a different conditional population. Compare pairs WITHIN it; do
not infer gain/attenuation from absolute scores versus primary or old nhanes04.
LLOD/sqrt2 fills remain censored reported substitutions, not measured low values.
No winner/threshold/transform selection should follow a favorable sensitivity.

## Required pre-run/readout checks

1. Current revised SSCARD schema/full marginal source counts, corrected weight
   columns and matching historical source hashes; no small-subdomain count forced
   to equal an entire-release marginal frequency.
2. One-to-one SEQN and correct1999/2001/2003 cycle reconstruction; masked design
   support/allowed codes; full DEMO versus eligible PSU inventory gate.
3. W6 exactpiecewise rule, no per-cycle scoring normalization; finite positive
   eligible weights, same source weights in every primary arm.
4. Common-mask equality, detection flags, sentinel/zero decoding and complete
   source/endpoint gates. Actual source0 bits must be repaired exactly; no epsilon
   rounding of tiny values into flag0.
5. Core/assay source definitions pinned before outcomes; means/scales and lambda
   choices fitted inside every actual nested training split and bootstrap.
6. Interpretation: this is pooled reused-calendar internal assessment. Storage,
   assay and baseline-era effects may remain; random whole-PSU CV across pooled
   cycles is not forward temporal or external validation. No new cycle-indicator
   variant or partition search is requested by this review.

No repository file was edited by this reviewer. The source loader still requires
its own implementation/test audit before empirical source qualification is
complete; no individual-data or result verification is claimed by this memo.

## Addendum: source loader now reviewed before empirical execution

After the initial method review, read the completed cohort.py and retrieve.py.
This supersedes the earlier statement that the loader was absent. No real files
or comparison metrics were opened. SOURCE_SPEC and locked metadata remain
pre-empirical declarations, not a claim this reviewer retrieved raw bytes.

Source verification precedes table parsers and binds schema/plan/lock, all file
names, exact official URLs, byte counts and SHA. Every scientific DAT/XPT must
also match the prior nhanes04 fingerprint list; merely changing a self-reported
manifest cannot authorize a different historical release. SSCARD XPT comes from
one unique official catalogue anchor, with no URL-guess fallback. Missing/unsafe
or ambiguous anchors, failed files and source rebinding stop the import.

The loader checks corrected full-release schema and21890row count, complete
marginal detection/pristine/weight counts, aligned result/flag missingness,
positive finite reported results, one-to-one cross-cycle keys and mortality/demo
universes. Source numeric-type checks and the root's finite-positive-integer XPT
key guard precede joins; exact source-byte zero repair remains upstream.

The six-year rule uses W4 for cycles1/2 and W2 for3, independent of unrelated
missing W4 in2003. It preserves source weights and uses no arbitrary cyclewise
rescaling. Primary flag00 and sensitivity flags01 are separate shared masks;
flag2 cannot enter either. Raw released source replacements retain their stated
interpretation; the core does not re-impute undetected proteins.

Before model arrays are constructed, full DEMO and eligible masked-pair sets
must match exactly and every retained stratum must have>=2units, allowing3.
The root now exports both masked-design inventories as public design metadata,
not participant IDs. Counts alone are no longer the only record of this gate.
No arbitrary missing-unit/singleton workaround was found.

Implemented only tests/test_b2m31_sources.py:25 synthetic/mock checks passed.
Tests cover official unique links/URL safety, historical pins before parsing,
tampered self-consistent manifests, source count/schema/marginals and fake tiny
nonzero detection flags, piecewise pooling weights, flag00versus01 domains,
PSU3 preservation, dropped pair and singleton rejection. Fixtures are generated
software examples and contain no public cohort records or real model outputs.

No critical blocker found in this completed source implementation. Empirical
retrieval, actual hashes/counters, design inventory,8year support and test/CI
outcomes still must succeed; this review does not substitute for those checks.
