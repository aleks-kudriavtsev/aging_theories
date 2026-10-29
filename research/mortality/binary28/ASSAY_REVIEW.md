# Assay interpretability review for iteration28 (2026-10-29)

Metadata/primary-literature review before comparative model results. No individual
records, ZIP, accounts, new agreements or holder contact. No repository files
changed. This review supports a limited within-cohort **reported-values** analysis,
not analytical validation of the stored specimens or clinical concentration cutoffs.

## Primary sources and findings

### Original cystatin C: original measurement, not automatically modern calibrated Cys

[CDC original27a codebook](https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.pdf)
identifies SSCYSTAT as mg/L, first-exam surplus serum measured by Dade Behring
N Latex/PENIA on BNII. Published imprecision is assay-level, not per-record QC.
The prose range is0.23–7.25mg/dL, while the exact field label is mg/L and released
values span0.05–8.5mg/L. This is an unresolved internal metadata inconsistency;
no automatic x10, threshold exclusion or QC_PASS label is justified.

[Selvin et al., AJKD2013, DOI10.1053/j.ajkd.2012.09.013,
PMID23177702](https://pubmed.ncbi.nlm.nih.gov/23177702/),
[author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC3703771/), documents
NHANES cystatin calibration drift. Indexed primary full-text snippets give the
NHANESIII Deming equation Y=0.022+0.800*(CCRL2006 cystatin C), followed by factor1.12
for ERM471/IFCC traceability. Full table/body retrieval hit PMC reCAPTCHA in this
review; do not activate an exact calibration from this memo alone. Original
SSCYSTAT should remain explicitly labelled original reported values. This does
not prohibit fitting new coefficients within the same original-assay cohort.

### First-exam B2M/BTP: identity/units known, individual quantification flags unknown

[CDC36a release](https://wwwn.cdc.gov/nchs/data/nhanes3/36a/SSNH3BTP.htm) explicitly
labels SSB2M and SSBTP mg/L, Siemens N Latex immunonephelometric measurements
in2009 on stored original serum. The shared generic prose range again says
mg/dL, inconsistent with exact variables and QC means. Neither a repair factor
nor per-person LOD/LOQ/PASS/dilution flags is established. Do not turn numeric
presence into validated quantification or reject out-of-prose-range values as
proven assay failure. The same official page documents protein-subsample
weights; ordinary MEC weights are not a substitute.

[CDC37a repeat-exam release](https://wwwn.cdc.gov/Nchs/data/nhanes3/37a/SSNH3CYS.htm)
uses a different repeat visit and calibration context and has no population
weights. It must not be pooled with first-exam27a/36a or treated as paired
technical validation of every first-exam participant.

### Freeze/thaw evidence is reassuring but limited

[Juraschek et al., Clin Biochem2012, DOI10.1016/j.clinbiochem.2012.02.027,
PMID22425605](https://pubmed.ncbi.nlm.nih.gov/22425605/),
[author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC3575740/), compared
172 NHANESIII specimens before/after one additional freeze–thaw. Rank
correlations were0.90(BTP) and0.99(B2M); BTP decreased by mean0.05mg/L, B2M showed
no detected mean difference. Samples were stored at−70°C; remeasurement was
4–13months later. This supports persistence of a substantial reported signal,
particularly B2M. It does not identify each specimen's storage history, establish
absence of degradation over the original long storage, or validate clinical
accuracy at low/high concentrations. A mean study-level difference is not an
individual correction to apply to all records.

The independent author/institution record corroborates the primary abstract:
https://pure.johnshopkins.edu/en/publications/the-effects-of-freeze-thaw-on-%CE%B2-trace-protein-and-%CE%B22-microglobuli-4

## Decision for28

**Proceed only as limited internal prognostic analysis of released measurements.**
The three variables have sufficiently identified molecules, baseline specimens,
methods and published values to investigate whether their within-cohort measured
signal improves a separately named model. Per-record analytical usability remains
unverified. The analysis cannot validate a biological mechanism, kidney-function
formula, absolute clinical threshold or assay interchangeability, nor claim that
improvement necessarily reflects unconfounded biology instead of storage/batch
or selection effects.

No frozen1999+ coefficient should be transported onto these measurements. A new
NHANESIII core fit on original data avoids that transport operation; it does not
retroactively solve the assay/QC limitations.

## Prespecifiable transformations and sensitivities (methodological recommendation)

1. Primary: retain strictly positive finite reported values, log2 transformation;
   training-only centering/scaling and fixed training-only tail handling if used.
   No outcome-selected cutoff, no imputation of absent proteins. Exact same
   complete domain/folds for core and core+all3; primary contrast fixed before
   modelling. The claim is incremental prediction beyond the chosen core only.
2. Sensitivity with train-only empirical CDF/rank representation for proteins,
   scored on the identical participants/folds. It is invariant to common strictly
   monotone recalibration or unit scaling; use a fully prespecified tie rule and
   extrapolation at training extremes. It does not correct specimen-specific
   measurement error, batch-specific shifts or rank changes.
3. If primary log tails are not clipped, a second prespecified sensitivity using
   training-only winsorization at1st/99th percentiles can test influence of extremes.
   These are robustness quantiles, not LOD/LOQ or reference intervals. Never
   select the favorable transform and rebrand it the primary analysis.
4. Optional study-level BTP +0.05mg/L perturbation may be labelled a freeze–thaw
   stress test, not data correction. It is not necessary if two sensitivities
   would overburden the comparison; rank/tail handling already target main
   identifiable scale/extreme-value concerns. Avoid many post-hoc variants.

Every transform, scaling and penalty search must be fitted inside each inner
training fold, repeated in outer training fits and bootstrap refits. Using a
full-sample rank or full-sample clipping threshold would expose heldout data.
A log2 global unit multiplication becomes an additive offset removable by
training centering, but a nonlinear/affine calibration can change log spacing;
this motivates rank sensitivity rather than a claim of calibration immunity.

## Required reporting limits

State that numeric eligibility is not QC/LOQ eligibility; count missing/nonpositive
records explicitly. Keep all first-exam variables original, with sample weights
and time origin intact. Report sampling/complete-case selection, original assay
calibration and storage uncertainty. Bootstrap uncertainty conditional on reported
measurements does not include assay calibration or biological measurement error.
Clinical validity and independent external validation remain unestablished.
