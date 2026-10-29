# Cause33B — read-only pre-result code review

Дата: 2026-10-29. Protocol commit:
`d772b874198a9f7739225bbce5d1b6b6dfcd0c08`.
Проверены четыре новых scientificPython файла, соответствующие contract/plan,
новые тесты иworkflow. Person-level sources, новые predictive metrics,
HTTP/GET, empirical fits/bootstrap не читались/не запускались. Репозиторий
не изменялся. Число scientific modules в дереве —103 (99старых+4новых).

## Заключение

**Критических pre-comparison противоречий или блокеров в рассмотренном коде
не найдено.** Научная корректность orchestration соответствует зафиксированному
протоколу: отдельные002/001 labels, тот же общий domain/W6/PSU folds, один draw
на все восемь cause×arm fits, полное nested refitting и раздельная undefinedness.
Это code review до empirical run, не подтверждение будущих результатов/CI или
независимое воспроизведение participant-level анализа.

## Зафиксированные рассмотренные bytes

| Новый файл `research/mortality/cause33b/` | SHA256 |
|---|---|
|`__init__.py`|`01aa1d56485a9cc7805b95d05ad0a484dcce2d8b3f8b50a03d5b29d2fd25d9a0`|
|`retrieve.py`|`2b08cf53adfa7eb19052fab8c31d4e11f538027dd19aa66badef06e0978032b8`|
|`cohort.py`|`a723340edaf441cb2eeb02a0d6ffe24bbeb30a548c809bf3fce7a1f0c65af673`|
|`analyze.py`|`79d82765db720afd774a45c93750fe3e939282f88487154addd899e104d96f17`|

## Intake, source identity и численная projection

- `retrieve.check_study` проверяет pinned protocol/metadata и все99предыдущих
  scientificPython hashes. Wrapper привязывает plan, lock, exactmanifest bytes
  иfile→sourceSHA map до parser. Сохраняется прежний31intake, включая два
  `denied_not_retried` mortalityPDF; новых mirror/retry/ZIP путей не добавлено.
- `cohort.load` проверяет wrapper до `_load_prior`, сравнивает manifest после
  unchanged31loader, заново выполняет33A gate и сверяет всю scientific
  projection. Затем повторно проверяет wrapper и frozen33A bytes до возврата.
- Projection требует точных sets ключей, типов, списков, counts, statuses,
  code/fold fingerprints. Ровно preplanned`rtol1e−12/atol1e−9` относится к
  finite floating aggregates; boolean/int не становятсяfloat. Rawdata pins,
  native-to-W6 formula и measurement/domain проверки остаются точными до
  этой aggregate roundoff сверки. Proof экспортирует оба scientific hashes
  и фактические float-difference diagnostics, не скрывая несовпавшие bytes.
- Только exactproceed`['002','001']` и successful inheritedsource/unknown/
  support gates разрешают зависимые сравнения. Изменённые домены/коды/ранги/
  фолды не становятся «почти тем же» исследованием через tolerance.

## Отдельные labels, общий domain и отсутствие all-cause подмены

- `cohort.cause_targets` строит002/001 непосредственно изvalidateddead/cause/
  t_exam и исходногоcommon order; `analyze.run` использует именно эту функцию.
  Старый `b2m31.analyze._domain_arrays` сall-causey не вызывается.
- Смерти при0/96месяцах включаются. Известные earlycompetingdeaths,
  смертипосле96 иsurvivors>=96 являютсяnoncases выбранногоtarget. Early
  unknowncause блокируетобаtargets; код дляalive, undeclared/fractionalcause,
  неверныеstatus/time иearlycensoring не исправляются/не удаляются.
- В`run` повторно проверены6044/6546 ипредзаданныеevent/controlcounts доfit.
  Features/proteins/weights берутся в одномcommon order, proteins толькоlog2
  positivefinite reportedCys/B2M. `_validate_arrays` проверяет совпадение
  eligible→commonweights/clusters/folds и наличие обоихклассов длякаждогоtarget.
- Проверка непересечения двухeventindicators — исходная mutuallyexclusive
  underlyingcause семантика, не интеграция predictions. Отдельныеprobabilities
  нигде не суммируются; joint/FG иbetween-causeeffectcontrasts не реализованы.

## Обучение и целостность folds

`analyze` вызывает сохранённый `b2m31.model.nested_crossfit` отдельно для
каждогоtarget/arm. Code path:

1. Outerassignment наeligiblePSUs черезunchanged`b2m31:outer`, затемcommonmask.
   Complete eligible/common inventory33A сохраняет исходный87-PSU набор.
2. В каждомouter training inner3 строятся по`b2m31:inner:{outer}`. Cluster не
   может пересекатьboundary; scoringlabels не подаются вtraining/tuning.
3. Actualtraining учитmeans/populationSD/constantflags; λ1/10/100 выбирается
   поinner-heldoutBrier своегоtarget сlargerλtie rule; freshouterfit послеtuning.

`_check_shared_preprocessing` сверяетoutcome-free parameters одинаковогоarm/
outerfold междуcauses. Финальная wordingточна: проверены **доступные matched
fits**, unavailablepairsпослеfailure не объявлены проверенными. Это comparison
preprocessing, не cause×markerinteractiontest. Likelihood normalization и
ridge primitives31 не модифицированы; sourceW6 не заменены native32weights.

## Один500-draw stream, fullrefits и zero-weight rows

В`analyze_arrays` RNG создаётся один разсseed20261105. Каждый изровно500
циклов вызывает`design_bootstrap_weights(eligible,rng)` **одинраз**, затем
`wr=full_weights[common]` используется для002/001 ивсехчетырёхarms. Full
eligible maskeddesign, включая3PSUвstratum1, обслуживается прежним m−1 /
m/(m−1) rescaledrule. Source/evaluationcaseweights не определяют draw.

Каждый`_fit_score` заново вызываетnestedcrossfit; не переиспользуются point
coefficients или pointpredictions. Исходные строки/clusters/folds сохранены,
включаяreplicatezero-weightrows. Новые draws не создаются дляreplacementfailures.
Вwideexport естьalignedweightbytes SHA256 и aggregateweight/eventweights,
позволяющиеаудировать pairing без раскрытия vectors.

## Failure semantics и интервалы

- Cause/armfiterror локален; вторичныйpair failure не стираетдоступныйprimary
  pair ипрогнозы другойcause. Fit/metricerrorreasons отражены отдельно.
- Preservedmetrics оставляетBrier приsingle-classundefinedAUC илокально
  undefined/nonconverged calibration. Calibrationfit не меняетscoredp.
- `flatten_scores` сохраняет оба pairedcontrasts; `paired_bootstrap.csv`
  содержит отдельныеBrier/AUC statuses, поэтому computedBrier не маскирует
  undefinedAUC. Все500attempts и500×2causes×2contrasts=2000pairedrows остаются.
- `intervals` считаетfinite-definedpercentiles сvalid/undefinedcounts.
  Дополнительныеflags явно говорятmarginalexploratory, notfamilywisecoverage
  инеучётendpointselection. Calibrationstatusreasons не маскируются нулём
  predictivefitfailures. Нетfixed-predictionCI вместоfullrefitCI.

## Export и фактическая граница review

Шестьoutputs —RESULTS_SUMMARY,FOLD_PARAMETERS,CALIBRATION_BINS,SOURCE_MANIFEST
JSON иwide/pairedCSV. Всеpayloads сериализуются до создания новогоoutputdir;
существующийoutput не перезаписывается. Pointouterfit records допустимы;
bootstrap coefficients, source records, participantIDs, individuallabels/
prediction/weightvectors не включены. Протокол разрешает толькоhashцелого
alignedweightvector какaggregatefingerprint. Markedstdoutblocks печатаютexact
сериализованныеfilebytes; hashesnonresultpayloads записаны вresult.

Pointpredictions остаются вRAM; point/statuses не объявляютобязательныйуспех
всехfits. `comparative_model_fitting_attempted=True` и
`comparative_models_fitted=bool(parameters)` отделены отinheritedprefitaudit.
При empiricalreporting всёравно необходимо смотреть конкретныеpoint/replicate
statuses, а не считать processsuccess научной гарантией.

Прочитаны3intake,17cohort/endpoint/projection и10analyzer testdefinitions.
Новыеanalyzerтесты —syntheticorchestration сmockedfit/metrics иrealpure
preprocessing; они не являютсяновой numericalvalidation самихridgefits.
Часть500-attempttest намеренномокируетfitting, что позволяет проверить
pairing/failurecount/export без новыхempiricalcomputations. Старые31numerical
primitives сохранены иихproductionbehavior проверитещёscientificempiricalCI.
Вэтомread-onlyreview тесты не выполнялись; root/агенты отдельно сообщают
3+17+10passingtests, полныйsuite/remoteCI должен бытьпроверенroot.

В scope нетlabreportv3/presentation, deployment, merge илиsourceaccountaccess.
Ни новыхpredictiveзначений, ниclaimsополезностиB2Mдлякакой-либопричины этим
code review не получено.
