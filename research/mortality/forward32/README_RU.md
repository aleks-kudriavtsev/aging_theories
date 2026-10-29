# Итерация 32: ретроспективная проверка B2M по календарю baseline

Строгое разделение по baseline: **обучение 1999–2002**, оценка только
**2003–2004**. Исходы поздних участников не входят в подбор, preprocessing
или обучение. Все периоды уже использованы программой; это exploratory
stress test, не untouched holdout, внешняя валидация или перенос frozen30.
Восьмилетние early labels становятся известны в 2007–2010, после baseline
late2003–2004: модель нельзя представлять как доступную для применения в 2003.

**Основной вывод ограниченный:** point Brier улучшается при добавлении B2M
к routine, но верхняя joint-refit граница всего−0.0000130815 от нуля.
При 500 Monte Carlo draws это не устойчивое подтверждение эффекта. AUC interval
пересекает 0. Убедительная добавочная ценность B2M сверх routine+Cys по основному
joint-refit интервалу не установлена. Узкий fixed-training интервал её не заменяет.

## Один заданный late domain и четыре фиксированные модели

Common early: **4020 участников /383 смерти** за 8 лет; common late:
**2024 участника /203 смерти**. Source/feature/flag00 критерии неизменны 31:
screening age40–75, известный exam-origin8year статус, полный routine11 и оба
положительных Cys/B2M с флагами 0. Внутри каждого периода один domain для всех 4
arms. Нет новых подгрупп, transforms, горизонтов или выбора победителя панели.

| Поздняя модель | Brier | AUC | Средний риск | Calibration slope |
|---|---:|---:|---:|---:|
| Routine | 0.05962233 | 0.80248093 | 7.639543% | 0.988717 |
| Routine + B2M | 0.05787822 | 0.81411187 | 7.657141% | 1.027761 |
| Routine + Cys | 0.05865295 | 0.81364399 | 7.066619% | 1.006222 |
| Routine + Cys + B2M | 0.05790775 | 0.81589992 | 7.323026% | 1.025284 |

Наблюдаемая взвешенная поздняя смертность 7.132987%. Уверенное улучшение
калибровки не заявляется. Mean risk и диагностические slopes не использовались
для исправления прогнозов; late intercept update/recalibration отсутствуют.

## Основная неопределённость: joint_refit

Первичный контраст сохранён 31: B2M сверх routine; вторичный — B2M сверх
routine+Cys. Для каждого ровно 500 парных попыток, без redraw. Mode distributions
не смешиваются и не выбираются по узости/знаку интервала.

| Контраст | Point ΔBrier | PRIMARY joint-refit95% percentile | Fixed-training95% percentile |
|---|---:|---|---|
| B2M сверх routine | -0.0017441124 | [-0.0032533156; -0.0000130815] | [-0.0031857193; -0.0002790308] |
| B2M сверх routine + Cys | -0.0007452084 | [-0.0015307437; 0.0000972323] | [-0.0013643449; -0.0001108726] |

| Контраст | Point ΔAUC | PRIMARY joint-refit95% percentile | Fixed-training95% percentile |
|---|---:|---|---|
| B2M сверх routine | 0.0116309388 | [-0.0027032630; 0.0256862840] | [-0.0021557249; 0.0237208418] |
| B2M сверх routine + Cys | 0.0022559315 | [-0.0062030549; 0.0109507353] | [-0.0057379289; 0.0087951465] |

Fixed-training условен по единственному fitted early model и отражает только
позднее PSU-resampling. Joint-refit независимо resamples early+late PSUs,
повторяет early preprocessing,3fold tuning и full early fit. Одинаковый late
draw используется в обоих режимах и всех arms; RNGlate20261103 и
RNGtrain20261104 независимы. Нельзя подменять неопределённость обучения
fixed-mode интервалом, особенно для secondary: fixed исключает 0, joint включает 0.

Все Brier/AUC500/500 определены, ошибок прогнозных refits0 в каждом mode/arm.
Все joint calibration diagnostics500/500. В fixed mode routine+Cys
intercept/slope499/500,1undefined: их defined-only bounds условны по этим fits,
не полное nominal95% покрытие. Остальные fixed diagnostics500/500.
Bootstrap условен по формам, source measurements, fixed folds и истории выбора;
не включает assay/storage/selection-history/calendar-future uncertainty и не
является официальным NCHS variance estimator. Граница near0 не оправдывает
повторный выбор seed, увеличение resamples ради знака или p-value discovery.

## Правильные веса и достаточность наблюдения

**TrainWTSSCB4Y**, **lateWTSSCB2Y**, непосредственно из исправленного
SSCARD_AJuly2022. ПолеW6 из 31 явно заменено, исходные веса/W6 оставлены только
для provenance. Метрики разных периодов не объединяются. Это отдельные
received complete-case/measurable domain estimands; веса не устраняют
дополнительный linkage/missing/flag отбор. Нет безусловной national validation.

| Показатель | Early1999–2002 | Late2003–2004 |
|---|---:|---:|
| Common n /events8 |4020 /383|2024 /203|
| Masked strata /PSU |28 /57|15 /30|
| Event-bearing PSU |57|30|
| Kish effective n |2572.544991|1333.692684|
| Native common weight sum |99016388.114702|101970999.131736|
| Native common /eligible coverage |93.754349%|91.400258%|
| Minimum survivor followup,months |201|177|
| Non-pristine /pristine |33 /3987|10 /2014|

Поздний WTSSCB4Y отсутствует у всех 2024common участников; это ожидаемо и не
удаляет их. Его sum(skipna)=0 не означает индивидуальные нулевые веса.
РаннийW4sum =1.5×старыйW6sum, lateW2sum =3×старыйW6sum. Полные eligible/common,
raw/squared/native/percycle sums экспортированы, coverage/Kish восстановимы.
Все full=eligible=common masked inventories сохранены; stratum1 имеет 3PSU,
остальные 2. Participant/PSU overlap между периодами 0. Это coded public masked
идентификаторы, не доказательство независимости истинной географии.

Gate>=20event rows и>=20nonevent rows с положительным весом в каждом периоде
пройден; это прагматический computational minimum, не power/EPV/calibration
certificate. При провале pipeline сохраняет blocked-feasibility без fits.
Eligible-derived early foldmap57PSU/3folds по 19PSU экспортирован; равенство
eligible/common inventories доказывает неизменность hash partition. Zero-weight
resample rows сохраняются; innerfold mapping не переназначается по исходам.

## Раннее обучение и assay ограничения

Все 4early final lambdas =10, выбраны только в 3early inner folds из 1/10/100,
при tie1e-12 — большая lambda. Weighted mean/population SD обучаются только
в фактическом early training, потом новый full early fit; intercept не penalized.
Новые artifacts MODEL_ARTIFACTS содержат 4fit/36inner diagnostics/3fold support.
InnerBrier и objective — training diagnostics, не дополнительная validation.

Core/units/creatinine calibration и логарифмы точно 31. Никаких новой Cys
standardization, NHANESIII CEP corrections, imputation, clipping или late scaling.
Серия serum assay2018–2020 не равна biological baseline. Публичный SSCARD не
даёт batch/run/individual-assay-date/storage-duration IDs; SSPRIS не является
batch. Поэтому calendar, storage и assay effects неразделимы по этим данным.
Экспортированы заранее заданные outcome-free mean/SD summaries; они не
используются для коррекции или выбора вариантов. Флаг 0 не является QC/LOQ pass.

## Код, источники и независимый аудит

Protocol commit377778c и implementation clarification0a27717 предшествуют новым
метрикам. Финальный scientific head **59b0626767c0169107652cc7b7d51a416c8e7668**,
CI run**40917437691**: empirical и 54 теста Python3.11/3.13 success;
workbench40917437642success. Реальный PRmerge checkout
**def39d99ab86bb6f7b022554bdb9f7cb077e3d27**,95scientific file hashes равны
head. Старыйbf2run40917121566 superseded после исправления no-fit empty-period
пути; его метрики не читались/не выбирались как основание результата. Поздний
docs-only head не считается новым scientific CI и не повторяет расчёты.

Все 50 источников и 25historical data SHA проверены. Full observed SSCARD marginals
теперь экспортированы для всех flag/weight/pristine категорий, вместе с native
weight arithmetic. Два прежних mortalityPDF403 остаются prior-denied/not-retried,
отдельно от 50retrieved inputs; ни обходов, ни новых PDF requests нет. Сохранены
current source+wrapper manifest bytes; rawdata не экспортируются.

Независимый численный аудит **203/203** пересчитал все 80 интервалов (40 в каждом режиме)
из 500+500wide и 1000paired CSV, проверил одинаковые late draws, постоянные early weights в fixed и изменяемые
в joint, артефакты и source hashes. Независимый source/domain аудит —
**380 проверок /0 отказов**. Это аудит агрегатов/кода/провенанса: исходные
XPT/ASCII не перерабатывались аудиторами, индивидуальные Brier/AUC и gradients
не восстанавливались из rawrows. Нет нового clinical/external claim.

```bash
python -m pip install numpy==2.3.5 pandas==2.2.3 scipy==1.17.0 scikit-learn==1.8.0
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_forward32*.py' -v
python -m research.mortality.forward32.retrieve --out /tmp/new_forward32_inputs
OPENBLAS_NUM_THREADS=1 python -m research.mortality.forward32.analyze --source /tmp/new_forward32_inputs --out /tmp/new_forward32_results
python tools/verify_forward32_aggregates.py --repo . --log /tmp/saved_forward32_job.log --ci-jobs /tmp/forward32_independent_ci_jobs.json --ci-log-dir /tmp --out /tmp/new_forward32_audit
python tools/audit_forward32_sources.py --repo . --log /tmp/saved_forward32_job.log --reference research/mortality/forward32/SOURCE_EXPECTED_REFERENCE.json --out /tmp/new_forward32_source_audit
```

Для численного аудита нужны completed-job plaintext logs с именами
forward32_independent_job_JOBID.log. Все analysis out-каталоги новые. Точные
7JSON/CSV stdout blocks,source manifests,coefficients,bootstrap deltas сохранены
в repo; ZIP не скачивался, артефакт не единственное хранилище результата.
Точныеrun/job/SHA/expiration — RESULT_PROVENANCE.json.

## Актуальный scope после 32

Текущая all-cause оценка 32 сохранена. Следующая отдельно обсуждённая работа —
вклад маркеров по доступным определённым причинам смерти: ранжирование по
взвешенному mortality burden в явно заданной когорте/горизонте вместе с числом
событий; начальный ориентир до 5 поддержанных причин, без принудительного
включения heterogeneous other. Не рейтинг по AUC и не перенос ранга NHANES
на Россию/мир. Сначала code/year/weight/event gate, затем фиксированные
baseline+marker comparisons с корректным учётом competing outcomes.

Совместная интеграция новых причин и Fine–Gray **отложены**. В 32 их реализации,
benchmark или новых cause fits нет; был только read-only обзор, остановленный
на временных заметках. Исторический five-cause engine/прочие модели/main/UI
сохранены. Маркеры не назначаются исключительно одной причине, независимые
binary probabilities не суммируются. Новых scheduled tasks, accounts,
agreements,restricted requests,контактов,merge/approvals/deployment нет.
