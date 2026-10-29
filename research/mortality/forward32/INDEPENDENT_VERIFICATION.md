# Независимый численный аудит32

**203/203 PASS; общий исход: PASS.**

Научный commit `59b0626767c0169107652cc7b7d51a416c8e7668`; final CI run40917437691, empirical job123601312507. Superseded run40917121566 и его метрики не использованы.

Early n=4020,events8=383; late n=2024,events8=203. Event-bearing PSU:early57,late30.

Все80 интервалов пересчитаны независимо из500 строк каждого режима.1000строк paired CSV согласованы с режимами; late draws идентичны fixed/joint, early weights постоянны только в fixed.

| Контраст | ΔBrier point | Joint refit95% percentile | Fixed training95% percentile |
|---|---:|---|---|
| primary | -0.001744112378064315 | [-0.0032533155852408897,-1.3081464121040813e-05],valid=500 | [-0.003185719275323107,-0.00027903081763082545],valid=500 |
| secondary | -0.000745208357597095 | [-0.0015307437334940571,9.723232595999206e-05],valid=500 | [-0.0013643449011064806,-0.0001108725807720666],valid=500 |

Primary joint-refit верхняя граница **−0.000013081464** очень близка к0. При500 Monte Carlo draws это ограниченный исследовательский результат, а не устойчивое подтверждение. Для secondary отрицательный fixed-training CI не заменяет пересекающий0 основной joint-refit CI.

| Контраст | ΔAUC point | Joint refit95% percentile | Fixed training95% percentile |
|---|---:|---|---|
| primary | 0.011630938823259718 | [-0.002703263025251312,0.02568628403157289],valid=500 | [-0.0021557248918267656,0.023720841820909336],valid=500 |
| secondary | 0.002255931465422689 | [-0.006203054865342975,0.010950735286524446],valid=500 | [-0.005737928873712161,0.00879514653908126],valid=500 |

Metric-local calibration valid counts из500:

| Mode | Arm | Offset | Intercept | Slope |
|---|---|---:|---:|---:|
| fixed_training | routine | 500 | 500 | 500 |
| fixed_training | routine_b2m | 500 | 500 | 500 |
| fixed_training | routine_cys | 500 | 499 | 499 |
| fixed_training | routine_cys_b2m | 500 | 500 | 500 |
| joint_refit | routine | 500 | 500 | 500 |
| joint_refit | routine_b2m | 500 | 500 | 500 |
| joint_refit | routine_cys | 500 | 500 | 500 |
| joint_refit | routine_cys_b2m | 500 | 500 | 500 |

Ограничения:

- Retrospective baseline-calendar split, not operational prospective validation: early1999–2002 eight-year labels can mature in2007–2010, after late2003–2004 recruitment. Late labels are absent from fit, but training labels were not all available at the late baseline calendar.
- Both periods and late outcomes were previously used. This is exploratory reused-cohort evaluation, not untouched or independent external validation; no frozen30 transport or clinical readiness claim.
- Primary uncertainty is joint_refit. fixed_training conditions on one fitted early model and excludes training uncertainty; the two500-draw distributions are separate, not pooled or interchangeable.
- The primary joint-refit upper Brier bound is only about1.3e-5 below zero. With500 Monte Carlo draws and this near-zero boundary, describe the result as limited exploratory support, not robust confirmation or a formal discovery; do not replace it with the narrower fixed-trained interval.
- All80 percentile intervals are independently recomputed by sorted-order linear interpolation; finite-sample nominal coverage is not established for nonlinear penalized fitting, discrete tuning and AUC/calibration. Undefined metric-local draws remain disclosed.
- Bootstrap holds the feature forms, domain, inner partition and observed measurements fixed. Assay/storage/calibration/prior-selection and unknown future-calendar uncertainty are excluded.
- Individual assay batch/date/storage duration and release-specific Cys traceability are unavailable. Calendar differences cannot identify biological drift separately from storage/batch effects.
- Source bytes were fetched/checked by the empirical runner. This audit reconciles exported fingerprints, locked historical hashes, exact stdout blocks and source-wrapper hashes; it is not an independent second download or local rehash of raw public records.
- Only aggregates/source metadata/coefficients were read; no person records, XPT/ASCII/ZIP downloads, predictions or model refits. Point Brier/AUC and optimizer gradients cannot be independently recomputed from aggregates; checks cover arithmetic, source hashes, fitted layout, reported diagnostics and bin moments.
