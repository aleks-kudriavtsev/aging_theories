# Независимый численный аудит итерации31

Проверки: **221/221**; общий исход: **PASS**.

Научный commit: `f27c7b47157db7fbaf31dc998bbb5c55c1c4f7e8`; CI run40908831254, job123570909752. Независимый запрос GitHub подтвердил empirical success и66 synthetic tests на каждой Python3.11/3.13. CI checkout — merge ref с этим commit; научные хеши совпадают с опубликованными результатами.
Общий домен: n=6044, событий8лет=586, effective n=3906.197150; eligible n=6546; 43strata/87PSU.

Аудитор использует только стандартную библиотеку и сохранённые агрегаты. Проверены точные marked-blocks, в том числе по независимо полученному completed-job log,50 источников и25 исторических data SHA,12 metadata pins,84 сохранённых и6 новых научных модулей,500 строк каждого CSV, все40 интервалов, парная арифметика, параметры фолдов и калибровочные корзины.

| Контраст | ΔBrier | percentile95% CI | ΔAUC | percentile95% CI |
|---|---:|---|---:|---|
| primary | -0.0014094962 | [-0.002214549546082801, -0.00065692033150264] valid=500 | 0.0115566625 | [0.004359476262077203, 0.023543881102270993] valid=500 |
| secondary | -0.0004316326 | [-0.0010227331995214138, 0.00015128361195393316] valid=500 | 0.0022732590 | [-0.0027857527283261935, 0.00816266733786681] valid=500 |

Fold records: 40; sensitivity: computed, n=6054. Model failures: `{'routine': 0, 'routine_b2m': 0, 'routine_cys': 0, 'routine_cys_b2m': 0}`.

У калибровки меньше определённых bootstrap-реплик, чем у Brier/AUC:

| Модель | Offset valid/500 | Intercept valid/500 | Slope valid/500 |
|---|---:|---:|---:|
| routine | 498 | 481 | 481 |
| routine_b2m | 498 | 489 | 489 |
| routine_cys | 500 | 492 | 492 |
| routine_cys_b2m | 497 | 487 | 487 |

Эти калибровочные интервалы условны по определённым диагностическим fits. Нулевые ошибки прогнозных refits не означают500 определённых slope/intercept.

Ограничения:

- No raw participant/XPT/ASCII/ZIP downloads, individual prediction access or refitting occurred in this audit.
- Source SHA values are reconciled across locked historical pins, current exported manifest, runtime audit and success log; this is not an independent local rehash of raw source bytes.
- Full-release SSCARD marginal counts and corrected weight counts are checked by the pinned runtime guard; no separate observed full-release marginal table is exported. The independent audit verifies those guards and exported conditional-domain marginals, not a second XPT parse.
- Person-level Brier/AUC and fitted calibration regressions cannot be recomputed from aggregate exports; independently checked quantities are interval arithmetic, pair differences, metric identities and calibration-bin weighted moments.
- 500 rescaled-stratified-PSU bootstrap draws repeat training, scaling and inner tuning but hold the original whole-PSU partitions fixed. Assay standardization and prior model/marker selection uncertainty are excluded.
- Cohort1999-2004 and Cys have prior use; B2M was previously loaded. These are exploratory new binary models on reused persons, not independent external validation.
- Both official mortality PDFs remained denied HTTP403, denied_not_retried, separate from50 required sources. The amendment retains prior primary definition review, pinned public ASCII hashes and two-reader agreement without changing the endpoint.
- The flags0/1 sensitivity uses a different conditional domain and supplies within-domain paired point comparisons only; it is not an additional bootstrap interval or external validation.
