# Scientific validation 08

A participant-disjoint temporal test of the frozen historical-US all-cause model.

- Analysis was locked in commit f61b86c6fc00a43ae2490cf001b06c680271af63 before retrieving the new outcomes; see ANALYSIS_LOCK.md.
- Development: NHANES1999-2002, n4549. Calibration: previously examined 2003-2004, n2381. New evaluation: 2005-2008, n6049, 958 deaths within10years.
- Primary analysis includes frozen training imputation. Complete-profile sensitivity (n5287,760 deaths) was specified in advance. Missing-profile decomposition is post hoc.
- Full results, including increased Brier in the primary imputed sample, are in ../ITERATION_08_RU.md. Do not generalize the favourable complete-profile result to missing profiles.
- recalibration_spec.json is evaluated research configuration, not a new clinical model. The Workbench0.7 default and its frozen coefficients are unchanged.
- Use source checkout for the optional ../replay_temporal08.py module. The existing0.7 wheel does not package this new research specification.
- The accompanying mortality_science08 archive contains complete aggregate metrics, source manifests and scientific test logs. Person-level records are not redistributed. In git: scientific code, report, abstract, compact results and calibration specification.

Reproduction (new output directories):

```bash
python -m research.mortality.validate_temporal08 --old-source OLD_CDC_INPUTS --new-source NEW_CDC_INPUTS --model research/mortality/reanalysis06/model_bundle.json --out results08
python -m research.mortality.sensitivity_temporal08 --source NEW_CDC_INPUTS --model research/mortality/reanalysis06/model_bundle.json --results results08 --out sensitivity08
python -m research.mortality.audit_reader_reference08 --reference R_ReadInProgramAllSurveys.R --old-source OLD_CDC_INPUTS --new-source NEW_CDC_INPUTS --out reader_audit.json
```

New input URLs are specified in .github/workflows/nhanes-validation08-inputs.yml; reference reader in validation08-reader-spec.yml. Neither workflow has a schedule or repository-write permissions. Do not assume an expired artifact is a current download. The original file fingerprints accompany the full scientific archive. Changes to an input version must be recorded, not called the same validation.

The 1000 paired rescaled PSU bootstrap intervals condition on frozen training and calibration estimates; they exclude fitting/model-selection uncertainty. Clinical use, RU/DE calibration, cause-specific predictions and20year extrapolation are not enabled.

Conference text is in ABSTRACT_IYGF2026_RU.md; the one-page DOCX is supplied separately. The conference submission itself has not been sent.
