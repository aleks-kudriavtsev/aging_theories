# Независимый численный аудит33B

**265/265 PASS; общий исход: PASS.**

Protocol `d772b874198a9f7739225bbce5d1b6b6dfcd0c08`; scientific commit `ab09ca0b605b44544f49d21d1161c61609ef92d1`; CIrun40928005361,job123638594625; actual merge `819b97537bdbf0fec55f50e4aa4a882fd3da67f5`.
Synthetic Python3.11/3.13: 30 tests each. Domain n=6044; events002=176,001=148; knowncontrols5868/5896.

Проверены80marginalintervals,500общихwide rows,2000pairedrows,source33Aprojection,99старыхмодулей,40pointfoldrecords при8успешныхfits иcrosscausepreprocessing invariance.

| Cause | Contrast | ΔBrier | Marginal95% percentile | ΔAUC | Marginal95% percentile |
|---|---|---:|---|---:|---|
| 002 | primary | -6.816549640809022e-05 | [-0.00023702270516629177,0.00016682720868200362], valid=500 | 0.011313851887893578 | [-0.009970811649810526,0.039263552656707695], valid=500 |
| 002 | secondary | -3.869699335443816e-06 | [-0.00014822128519791506,0.0002113158154427012], valid=500 | 0.001098199291227897 | [-0.016781308558736267,0.016097964247292523], valid=500 |
| 001 | primary | -0.0002108011438316272 | [-0.0005336794627724806,0.00012534873092178572], valid=500 | 0.009343008260755847 | [-0.009427794829103017,0.04192884066555767], valid=500 |
| 001 | secondary | -9.088044701753561e-06 | [-0.0003101009316588021,0.00017543682314518846], valid=500 | 0.0031912647302395625 | [-0.007643984841639329,0.02526043725778365], valid=500 |

Это не familywise confirmation и не interaction test; разный знак/определённость по причинам не доказывает различия эффектов.

| Cause | Arm | Offset valid/500 | Intercept valid/500 | Slope valid/500 |
|---|---|---:|---:|---:|
| 002 | routine | 500 | 498 | 498 |
| 002 | routine_b2m | 500 | 497 | 497 |
| 002 | routine_cys | 500 | 495 | 495 |
| 002 | routine_cys_b2m | 500 | 499 | 499 |
| 001 | routine | 500 | 496 | 496 |
| 001 | routine_b2m | 500 | 497 | 497 |
| 001 | routine_cys | 500 | 499 | 499 |
| 001 | routine_cys_b2m | 500 | 498 | 498 |

Ограничения:

- Exploratory, internally nested comparisons on reused persons/proteins after full-cohort outcome-based source/support selection; no independent external or prospective validation.
- All80 intervals are marginal conditional percentile summaries, without familywise coverage, multiplicity-adjusted confirmation, p-values or a between-cause interaction contrast. Different interval signs or sizes do not establish different marker effects between causes.
- Controls include known early competing deaths, deathsafter96months and observed-through96 survivors. A competing death is a noncase of the fixed-horizon target, not censoring or proof of being alive.
- One eligible-frame PSU draw is reused across causes/arms; weight SHA values and moment/observed-risk invariants are checked, but the unpublished weight vector cannot be independently reconstructed or rehashed from aggregates.
- 500 refits repeat preprocessing/tuning/training with fixed rows/folds/forms. Metric-local undefined draws are retained; defined-only bounds are not complete nominal95% coverage. Assay/storage/calibration, public perturbation, endpoint/model-selection history and transport uncertainty are excluded.
- Source hashes, source projection, point parameter layout/diagnostics, bin moments, deltas and all interval percentiles are independently reconciled; no raw person-record read or model refit occurred. Person-level Brier/AUC and optimizer gradients are not recomputed from aggregates.
- Independent binary cause probabilities are not summed or advertised as a coherent joint model. Fine–Gray/joint integration, frozen models and laboratory reportv3/presentation remain outside this change.
