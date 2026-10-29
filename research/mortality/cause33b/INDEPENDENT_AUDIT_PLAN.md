# Cause33B: независимый аудит, подготовленный до предиктивных результатов

Дата: 2026-10-29. Рабочий план аудитора `/root/renal29`, вне репозитория. На этом этапе не читались новые33B прогнозы/метрики, не выполнялись raw/source GET, ZIP, обучение или Fine–Gray. Labreportv3 и презентация не затрагиваются. Protocol commit **d772b874198a9f7739225bbce5d1b6b6dfcd0c08** зафиксирован до прогнозов; LOCK20metadata/99historical modules иCONTRACT прочитаны. Final output schema/CI добавляются после готовности analyzer.

## Зафиксированный вопрос

Только две прошедшие33A причины: **002 — злокачественные новообразования; 001 — болезни сердца**. Домeн6044, горизонт96месяцев inclusive отMEC, W6=(2/3)WTSSCB4Y для1999/2001 и(1/3)WTSSCB2Y для2003. Это те же ранее использованные участники и прежние reported Cys/B2M, не новая независимая валидация.

Четыре прежних arms внутри каждой причины: routine, routine_b2m, routine_cys, routine_cys_b2m. Primary: B2M сверх routine; secondary: B2M сверх routine+Cys. Новая fixed-horizon binary cumulative-incidence модель не является cause-specific hazard, Cox, Fine–Gray или переносом frozen30.

Из frozen33A известно586 all-cause8y событий, cancer176, heart148, unknown early0. Для cancer known noncases=5868, из них410 competing early deaths; для heart5896, из них438 competing early deaths. Прочие5458 участников не имели смерти к96месяцам. Эти числа — арифметика уже сохранённого gate, не новые33B predictive results.

**AUC controls и Brier population включают все известные noncases**, в том числе competing early deaths. Competing death не считается censoring или доказательством выживания. Unknown early cause нельзя сделать negative: если gate/pins изменились, зависимый analysis останавливается. Поздняя смерть/подтверждённое выживание96месяцев — известный noncase этого горизонта.

## Нужные экспорты и проверки

1. Source/LOCK provenance: exact33A GATE_RESULTS/SOURCE_MANIFEST pins, ordered causes002/001, source50/data25, old99 scientific modules preserved, new complete code mapping,2 прежних PDF403 denied_not_retried. Источник/единицы/QC/LOQ ограничения прежние, без нового assay bridge или imputation. Final completed CI plaintext содержит точные JSON/CSV blocks; аудит сравнивает bytes/SHA с repo, не скачивает ZIP/raw.

2. Point metrics всех8 fitted models: n/events/noncases/competing-control count, weight_sum иweighted observed incidence, Brier/AUC/predicted risk/bias/OE/calibration diagnostics иmetric-local statuses. Одинаковые n/weights всех4arms иобеихcauses; observed incidence для каждой причины равна frozen33A weighted_events8/common_weight_sum. Два outcome-вектора disjoint, однако их модели не считаются совместной когерентной системой иих вероятности не суммируются.

3. **40 point parameter records**:2causes×4arms×5outer. Ordered features11/12/12/13, means/scales/constants, coefficients/intercept,3inner-Brier scores forlambda1/10/100, independent tie rule1e-12→larger lambda, gradient/convergence/normalized training weight. Same31 outer5/inner3 whole-PSU mappings: salt `b2m31:outer`; `b2m31:inner:{outer}`. Map задан наeligible inventory передcommon-mask;87PSU/43strata, stratum1 с3PSU. Не искать благоприятные folds.

4. Равенство conditional frames можно проверить сильнее одного флага: **training means/scales/constants для samearm/outerfold должны совпадать между002/001**, потому что preprocessing не зависит отy. Training/scoring n/weights/PSU тоже одинаковы. Event counts иcoefficients/lambda могут отличаться. Все preprocessing/tuning происходят внутри фактического training; точность оптимизации проверяется по diagnostics, без refits.

5. По зафиксированному протоколу: wide bootstrap CSV ровно500строк, оба causes,4arms,9metrics,arm status/error/calibration status,primary/secondary deltas. Парный narrow CSV —**2000строк `(cause,contrast,replicate)`**,500×2causes×2contrasts, один контраст настроку. Не удалять failed draw целиком из-за другого arm/cause и не redrawing. Индексы и common draw идентичны во всех8fits; seed **20261105**.

6. Для проверки shared draws экспортировать **общие** weight_sum, squared_sum, positive-weight rows/PSUs,eligible sums истабильный hash bootstrap weight vector без IDs/самого вектора; дополнительно по каждойcause positive-weight event/noncase counts иevent_weight_sum,allcause event weight sum. Во всехmodel metrics одногоdraw observed risk должен быть event_weight_sum/общийweight_sum. Позволяет проверить один и тот же resample, ане только assertion flag.

7. Ровно500 paired rescaled-stratified-PSU draws на full eligible frame, потом fixed common-mask. В каждойstratum выбиратьm−1PSU сreplacement иfactor m/(m−1), включаяm=3. Все actual training preprocessing,inner tuning иouterfits повторяются. Penalty-selection counts по8models суммируются до5×числоуспешных modelrefits; evaluation-only failure не равна training-fit failure. Fixed outer/inner identity сохраняется уzero-weight rows.

8. **80 интервалов**:2causes×(4arms×9metrics +2contrasts×2metrics). Независимый stdlib audit сортирует конечные значения иинтерполирует позиции `(n−1)*.025/.975`; проверяет lower/upper,valid/undefined,scope иdefined_draws_only. У calibration/AUC может быть меньше500 определённыхreplicates при успешных predictive fits. Conditional bounds поdefined draws не называть полным nominal95% interval.

9. Paired deltas каждой реплики пересчитываются expanded−baseline,primary иsecondary отдельно ивнутри одной причины. Проверяется narrow↔wide↔JSON interval correspondence. ΔBrier/ΔAUC междуразными causes не объявляются interaction; разный знак/CI significance сам по себе не доказывает difference междуcauses.

10. Calibration bins каждойpoint модели: суммыn/events/weights, weighted observed/predicted moments иbias identities восстанавливаются без персональныхпрогнозов. Персональные Brier/AUC иoptimizer gradients не могут быть независимо восстановлены изэтихaggregates: audit честно раскрывает данный предел.

11. Source33A reproduction projection: fields/типы/counts/folds точны, толькоfloat reduction roundoff допускает rtol1e-12/atol1e-9. Независимый recursive comparator неимпортируетcohort.py ипересчитываетfloating_values_checked/different/max_abs_difference иprojection hashes поэкспортированному scientific projection. Rawsource hashes иформулаnative→W6 остаютсястрогими; tolerance не применяется кцелымчислам,labels,выборке илифолдам.

## Интерпретация

Все95% интервалы **marginal exploratory**, без familywise coverage или multiplicity-free discovery. Есть2causes,2contrasts инесколькоdiagnostics на уже выбранных по событиям/историипредыдущихисследований моделях. Не делать вывод «маркер действует только дляоднойпричины», если CI другой пересекает0; не заявлять formal cross-cause interaction безотдельного заранее заданного контраста.

Абсолютные Brier/AUC разных причин имеют разную prevalence иcontrolcomposition, поэтому неранжируют причины по marker benefit. Приоритетпричин по burden уже зафиксирован33A; новыеpredictive metrics его неменяют. Независимые binary probabilities не суммируются, biological/causal/clinical interpretation ограничена.

Refit bootstrap учитываетусловную sampling/training/tuning variability прификсированных формах,folds,common domain иreported measurements. Неучитываетнеизвестнуюassay calibration/storage/batch, public mortality perturbation, исходнуюmarker/model selection, будущиеэпохи иливнешнийtransport. Gate support не является power/EPV/calibration certificate. Joint integration/Fine–Gray остаютсяотложенными.

## Порядок финального аудита

Послеimmutable protocol/schema → синтетические проверкиисследовательскогоpipeline → толькоfinal completed empiricalCI log → independent exact2JSON/CSVcapture → bytes/hashes/shared-draw/model/interval arithmetic → отчёт `/tmp/cause33B_audit.json/md`. Никакихповторныхfits/rawdownloads, новыхпорогов, seedsearch илиизменениястарыхрезультатов.
