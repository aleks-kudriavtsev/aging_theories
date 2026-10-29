# Independent routine29 aggregate and percentile audit

181/181 checks passed; 80 scientific Python files checked.

Common n=3730, events8=563; primary independent percentiles: {"delta_brier": [-0.0036370856653121127, -0.0008804441629611157], "delta_auc": [-0.004827672587886586, 0.031074852628574887]}.

- No participant records, raw source archives or individual predictions were accessed, and no model was rerun. Point Brier/AUC cannot be independently re-estimated from these aggregate exports.
- All 500 paired contrast rows were parsed independently and percentile endpoints recomputed by sorting and linear interpolation with Python standard library only.
- Calibration intercept/slope intervals for both primary arms have only 498 defined draws and 2 undefined; their endpoints were matched to stdout but cannot be independently recomputed from the primary-only CSV. They are conditional defined-draw intervals, not nominal full-500 coverage.
- Source fingerprints, SAS hashes and reader aggregate counts were matched to the successful job manifest; source data bytes were not downloaded locally. Additional direct official lab.sas metadata GET was denied by proxy (Tunnel connection failed: 403 Forbidden), so no independent local source-byte verification is claimed.
- Exported fit parameters cover main, partition and strict point refits. Per-bootstrap fit parameters are not exported; model convergence guards were reviewed and exported failure/penalty counts checked.
- Primary intervals condition on the fixed original PSU partition, complete profiles, reported assay values and specified model forms. Five additional partitions are exploratory sensitivities and are not pooled into those intervals.
- NHANES III is reused; this is internal research, with remaining historical assay/QC/LOQ and creatinine calibration uncertainty. Neither independent external validation nor clinical readiness is established.
