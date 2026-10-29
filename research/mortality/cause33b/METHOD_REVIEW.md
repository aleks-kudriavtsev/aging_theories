# Cause33B — методический review до прогнозов и comparative metrics

Дата: 2026-10-29. Scope: read-only чтение кода31/33A и документов. Новые
person-level данные, predictive metrics и HTTP/GET не использовались; fits,
bootstrap, сравнения и изменения репозитория не выполнялись. Совместная
cause integration/Fine–Gray, labreportv3 и presentation не входят в этот этап.

## Решение по дизайну

Предложенный отдельный fixed-horizon binary анализ методически допустим при
повторном прохождении точного source/cause/no-censor gate33A. Два заранее
выбранных исхода: злокачественные новообразования002 и болезни сердца001.
Переданные parent числа176/148событий — уже известные результаты support gate,
здесь они не пересчитывались. Они не являются гарантией power, точной
калибровки или определения всех bootstrap fits.

Общий домен остаётся6044участника, screening age40–75, baseline1999–2004,
96месяцев от exam, pooledW6,11routine terms и положительные reported Cys/B2M
с flags00. Это новая exploratory cause-specific задача на ранее использованной
выборке. Её pooled whole-PSU CV не является forward-calendar проверкой32,
untouched holdout, независимой внешней проверкой или переносомfrozen30.

## Что именно оценивается

Для каждой причиныk отдельно:

`Y_k = I(MORTSTAT=1 AND PERMTH_EXM<=96 AND UCOD_LEADING=k)`.

`E[Y_k|X] = P(T<=96,J=k|X)` — восьмилетняя cumulative incidence опубликованной
группы при полностью известных labels. Прямой binary likelihood допустим,
поскольку target фиксирован, а не предполагается неидентифицируемое «время до
целевого события при устранении остальных причин». Это не cause-specific
hazard ratio и не subdistribution-hazard модель. Раздельныеp_k нельзя
складывать: их сумма не обязана быть<=1, а1−sum(p_k) не является survival.
All-cause baselines31/32 и их показатели остаются отдельными и неизменными.

### Cases и единое определение controls для Brier/AUC/calibration

| Исход записи | ЗначениеY_k |
|---|---|
|Известная смерть отk в0…96месяцев включительно|1|
|Известная смерть от другой причины до96|0, competing noncase|
|Смерть после96, включая смерть отk|0 на выбранном горизонте|
|Assumed alive с follow-up>=96|0|
|Смерть<=96 с неизвестной причиной|Label не определён; gate блокирует анализ|
|Assumed alive с follow-up<96|Label не определён; gate блокирует анализ|

Именно ВСЕ известные noncases составляют AUC controls. Survivors-only controls
дали бы другой discrimination estimand. Конкурирующие смерти не удаляются,
не цензурируются и не помечаются как живые. При подтверждённом отсутствии
раннего цензурирования IPCW не нужен для этого binary target. Смерти при0 и96
месяцах не сдвигаются и не теряются. Unknown cause после96 не делает
восьмилетний label неизвестным.

Public2019 содержит perturbation follow-up/cause; оценка относится к
опубликованному endpoint. Истинная индивидуальная причина и error процесса
linkage/perturbation не восстанавливаются. Группа001 — болезни сердца, не
всяCVD и не отдельнаяИБС;002 — все злокачественные опухоли, не локализация.
MCOD-флагиDIABETES/HYPERTEN не являются underlying causes или predictors.

## Критическое место повторного использования кода

`b2m31.model.nested_crossfit` принимает произвольный валидный binaryy и
подходит как сохранённый numerical primitive. Но
`b2m31.analyze._domain_arrays` конструирует **all-cause**y. Новый wrapper обязан
построить и проверитьY_001/Y_002 явно; перенос этого старогоy без изменения
тихо воспроизвёл бы all-cause задачу. Нельзя использовать старые all-cause
коэффициенты как уже оценённую cause model.

Нужно повторно проверить source SHA, неизменённый общий mask, support/fold
fingerprints и отсутствие раннегоunknowncause ДО выражения`cause==k`:
обычное сравнениеNaN сk возвращаетFalse, что без этого guard ошибочно превращает
unknown death в noncase. Mask остаётся outcome-free и одинаков для двух причин
и всех arms. Никаких новых cause-based complete-case exclusions.

## Общие четыре arms и заранее направленные contrasts

Для каждой причины заново обучаются:

1. `routine`:11фиксированных routine terms;
2. `routine_b2m`:тот же core+log2reportedB2M;
3. `routine_cys`:тот же core+log2reportedCys;
4. `routine_cys_b2m`:тот же core+обаlog2proteins.

Первичный marker contrast в каждой причине —`routine_b2m − routine`;
secondary —`routine_cys_b2m − routine_cys`. Primary loss — weighted8yBrier;
отрицательнаяΔ означает улучшение. AUCΔ положительная при улучшении и
публикуется как сопровождающая оценка. Те же участники, признаки, исходные
веса и fold assignment во всех arms; no Cys-only/B2M-only winner selection.
Для reference normalization/risk skill нельзя выбирать удобный comparator
после результатов. Порог классификации и новые subgroup/transform поиски не
требуются этой задаче.

## Nested whole-PSU training и защита от утечки

Сохранить outer5folds по87eligible maskedPSUs, salt`b2m31:outer`, и inner3
внутри каждого outer training, salt`b2m31:inner:{outer}`. Страта1с3PSUs
остаётся в design; identity — `(SDMVSTRA,SDMVPSU)`, не один повторяющийся
PSU label. Конструкция по eligible design доcommon filtering и требование
совпадения common/eligible inventories33A обеспечивают те же mappings.

Каждая причина/arm/outer fit заново учит weighted mean/populationSD только
на actual training; constant feature scale1. Lambda1/10/100 выбирается по
pooled inner-heldout weightedBrier только своего target; ties1e−12 выбирают
большуюlambda. Затем preprocessing/coefficients обучаются на всём outer
training, и outer evaluation остаётся только для scoring. Признаки и penalties
не выбираются по outer target outcomes. Внутренняя calibration fit — только
диагностика готовыхOOF probabilities, никогда не изменение этих probabilities.

PooledW6: `(2/3)*WTSSCB4Y` для1999/2001, `(1/3)*WTSSCB2Y` для2003.
Отдельная normalization likelihood кtrainingn сохраняется как numerical
правило. Нельзя незаметно использоватьnative32W4/W2, class balancing или
другие веса случаев/controls: они меняют calibration/estimand.

## Ровно500 paired полных refit попыток

Рекомендуется один заранее записанныйseed и одна последовательность masked-PSU
множителей для обоихtargets и четырёхarms. На каждой попытке один и тот же
draw взвешивает все восемь fitted models; resampling выполняется по полному
eligible design до общегоmask. Дляstratum сmPSUs —m−1 выборок с возвращением,
factor`m/(m−1)` и сохранение исходных cluster/fold identities.

В каждой попытке для каждогоtarget/arm повторяются preprocessing, inner tuning,
outer full fits иOOF scoring. Фиксированные только predictions/coefficients
не дают training-uncertainty interval и не заменяют этот план. Нулевые replicate
weights не удаляютPSUs из fold mapping. Не повторять случайные draws до
успешного результата и не перетасовывать partitions ради convergence.

Армовая ошибка сохраняется локально дляcause/arm. Если secondary arm не
оценился, это не стирает определённый primary contrast другой пары или другой
причины. Отдельно учитываются model-fit failures, contrast failures и
metric-local undefinedAUC/calibration. Brier может быть определён, когда AUC
или calibration не определены. Все500 attempted rows сохраняются. Percentile
2.5/97.5% сvalid<500 — bounds условно на defined draws, а не полныеnominal95%CI.

Эти refit bounds учитывают declared resampling/training/tuning при фиксированных
folds, исходах/формах и выбранныхcauses. Они не покрывают ranking/support
selection, прежний выборB2M, assay/storage/calibration, source perturbation,
complete-case bias или внешнюю переносимость и не названы официальнымNCHS
variance estimator.

## Множественность и пределы интерпретации

До прогнозов назвать family: два primaryBrier contrasts (002/001), два
secondary beyondCys contrasts; AUC/calibration — сопровождающие outputs.
Публиковать все четыре, независимо от знака. Если используются обычные
per-comparison95%percentile bounds, явно сказать, что они не обеспечивают
family-wise coverage и анализ exploratory. Не выбирать благоприятную причину
или AUC вместо primaryBrier и не выдавать CI exclusion за confirmatory testing.

Причины выбраны после outcome burden/support gate, не по marker effect, однако
эта selection history всё равно является ограничением. Числа176/148 и
проходной gate не разрешают расширить список редкими причинами постфактум.

Каждый marker может участвовать в обеих моделях; коэффициенты — условные
ассоциации с вероятностью конкретного опубликованного исхода. Они не назначают
B2M илиCys «исключительно» одному заболеванию. Разность знаков/значимости
коэффициентов либо AUC/Brier между причинами не является cause×marker
interaction test. Разные penalties, распространённость исхода и case/control
состав дополнительно мешают такому сравнению. Эксклюзивность, causal biology
и interaction в этом дизайне не заявляются.

## Компактный checklist того, что сохранить до/после запуска

**До запуска:** immutablecomparisonplan+commit/lock; exactcausecodes002/001;
96месяцев и unknownpolicy; общийmask/weights/source/gate hashes; четыреarms,
два contrasts наcause, fixedfold fingerprints, lambda/tie rule;500draws/seed;
failure/undefined/multiplicity rules; новые output каталоги. Сравнительные
результаты не заменяют gate provenance и не перезаписывают31/32.

**Агрегатные diagnostics/results:**

- n,raw/weighted target events и известныеnoncases, observed weightedCIF,
  event-bearingPSUs/strata и effective event support;
- source/W6/mask/fold fingerprints, gateunknown/time checks и protocolcommit;
-40point outer-fit records:2targets×4arms×5folds, сλ,innerBrier grid,
  training/scoring support, means/SD, coefficients/intercept, convergence;
- pertarget/arm Brier,AUC,observed,predicted,biaspp,O/E,offset calibration,
  full intercept/slope/status и пятьweightedrisk bins сties policy;
- все500paired wide rows, model/contrast failures, valid/undefined counts,
  pairedΔBrier/ΔAUC и declared conditional intervals;
- scientific code/library hashes, независимая сверка агрегатов/CSV и CI scope.

Ни source records, ни participantIDs, ни labels/prediction vectors в exports
или stdout не нужны. Методологического terminal blocker при сохранении этих
условий не найдено; фактическая корректность новогоwrapper/targetкодировки и
refit pairing должна подтверждаться тестами до empirical execution.
