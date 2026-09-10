# Competing mortality 09

Historical US research model with coherent cause-specific cumulative incidence.
Development: complete profiles NHANES2005-2008 n5287; participant-disjoint test2009-2010 n3182.
No redistribution of participant data. Old Workbench defaults and old all-cause models are unchanged.

`ANALYSIS_LOCK.md` was committed before new input retrieval (`99e8c150`).
`model_bundle09.json` contains all trained preprocessing, coefficients and interval baselines.
The accompanying archive contains full `metrics.csv`, `intervals.csv`, `calibration_diagnostics.csv` and `calibration_bins.csv`; GitHub carries the compact `VALIDATION_SUMMARY.json` and full report.
The archive table `source_cause_counts.csv` preserves the ten original CDC groups separately from three model heads.
The archive `SOURCE_MANIFEST_*.json` identify official data versions; the repository retrieval workflow supplies URLs. An expired Actions artifact is not a permanent data source.

Three groups: diseases of heart (UCOD001), malignant neoplasms(002), all remaining/unknown deaths.
Do not call001 exclusively IHD or all cardiovascular disease. Do not call the residual natural/noninfectious.
Horizon5years primary;1/8years secondary. Full clinical AND routine profile required for both comparison panels.
1000 paired rescaled PSU bootstrap replicates are conditional on fitted models and omit training uncertainty.

```bash
python -m research.mortality.competing_risks09 --development-source DEV_INPUTS --test-source TEST_INPUTS --frozen-model research/mortality/reanalysis06/model_bundle.json --out new_results09
python -m research.mortality.diagnostics_competing09 --source TEST_INPUTS --model research/mortality/competing09/model_bundle09.json --out new_diagnostics09
python -S -m research.mortality.replay_competing09 --demo --out synthetic09.json
```

Source checkout required; existing0.7 wheel does not package this new research module/model.
Scientific computations need numpy/pandas/scipy/scikit-learn; tests additionally statsmodels.
Replay is standard-library only. All311 local tests used the original official inputs; CI unit tests alone are not that complete scientific run.

See `../ITERATION_09_RU.md`. New conference text is separate from the earlier08 text.
