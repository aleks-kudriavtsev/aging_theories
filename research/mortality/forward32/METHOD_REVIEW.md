# Итерация32: независимая проверка дизайна до новых сравнительных расчётов

Дата: **2026-10-29**. Reviewer: отдельный исследовательский агент `/root/renal29`.
Область проверки: предложение parent о strict forward train1999–2002 / late2003–2004, уже имеющиеся primary-source reviews итерации31 и начальный чистый `forward32/model.py`. Новые исходные данные, исходы/метрики32, live job logs и сетевые GET не читались; модели не обучались; файлы репозитория не изменялись.

**Заключение: условно методологически допустимо; непреодолимого blocker по известным определениям, весам или исходу нет.** Это новое строгое календарное обучение/оценивание на ранее использованной когорте. Оно может проверить устойчивость результата31 при запрете поздних данных в обучении, но не создаёт untouched holdout, независимую внешнюю валидацию либо проверенную клиническую assay-панель. До новых метрик необходимо зафиксировать приведённые ниже условия в протоколе/LOCK; это научная рекомендация, не approval PR или deployment.

## 1. Что именно оценивается

- Обучающая популяция: ранее определённая eligible NHANES1999–2002, screening age40–75, правильный положительный subsample weight, MEC/linked2019 mortality, затем зафиксированный полный routine11 + оба положительных reported proteins и оба detection flags0.
- Поздняя целевая популяция: аналогичный **условный общий домен** NHANES2003–2004. Веса не превращают его в всех американцев40–75: complete-case, measurable-protein и linkage ограничения сохраняются.
- Четыре неизменных arms: routine, routine_b2m, routine_cys, routine_cys_b2m. Один и тот же поздний домен/вектор весов для всех четырёх; в раннем обучении также общий домен для всех четырёх. Нельзя расширять comparator за счёт людей без белков.
- Primary estimand — late survey-weighted разность Brier8y `routine_b2m − routine` для **зафиксированной процедуры раннего обучения**. Secondary — `routine_cys_b2m − routine_cys`; AUC и калибровка остаются сопровождающими диагностическими показателями. Не выбирать лучший белок/arm/horizon по late метрикам.
- Strict chronology относится к биологическим baseline exam1999–2002 →2003–2004. Белки измерены retrospectively в2018–2020; это не проспективное внедрение теста в2003 году.
- 2003–2004 уже участвовал в прежней программе и31. Новый commit до расчёта32 обеспечивает порядок именно этой итерации, но не устраняет историю выбора маркера/формы по использованной когорте.

## 2. Веса и дизайн: корректный выбор и обязательные gates

Уже проверенные первичные документы:

- [NCHS Weighting Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx): исключение первых циклов1999–2002 и четырёхлетние веса; правило распространяется на subsamples.
- [SSCARD_A codebook, correctedJuly2022](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.htm): именно его пересозданные WTSSCB4Y/WTSSCB2Y.
- [NCHS Variance Estimation Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx): masked strata/PSU из DEMO, уникальная пара strata+PSU, сохранение design/domain логики.

**WTSSCB4Y для train1999–2002, WTSSCB2Y для test2003–2004 — обоснованный выбор.** Не объединять эти native веса в одну train+late метрику. Не нормализовать1999 и2001 отдельно: это изменит их относительную массу. Общая нормализация likelihood weights к n фактического раннего training split допустима при неизменной ridge convention; late scoring использует свой native survey weight.

W6 итерации31 отличался от соответствующего native веса лишь константой внутри каждого временного partition. Поэтому weighted means/Brier/AUC и training-normalized ridge должны быть инвариантны к такой общей константе; полезен синтетический тест этой инвариантности, а не новая подгонка на late.

До primary-mask: разделить полный qualified eligible frame по фиксированному SDDSRVYR/циклу, затем отдельно проверить **full-DEMO→eligible сохранность фактического pair inventory**. Early:28 strata, late:15 strata. В каждой eligible stratum≥2 PSU; первый ранний цикл допускает PSU3. Пределы кода не определяют число реально непустых пар — сверять runtime inventory. Полный design frame используется для bootstrap ДО common-case mask; не превращать исчезновение primary-case PSU в уменьшение m при rescaling. Identity — `(SDMVSTRA,SDMVPSU)`, не один повторяющийся SDMVPSU.

Эти masked units не являются доказательством независимых истинных географических PSU между эпохами. Независимые RNG-потоки — свойство алгоритма resampling, а не доказательство географической независимости выборок; возможная недоступная межпериодная зависимость остаётся пределом публичного дизайна.

## 3. Обучение и физическая защита от late leakage

Внутри early eligible/primary domain фиксировать hash-based whole-PSU3fold assignment без исходов. Один salt, без поиска лучшего разбиения. В каждом inner training отдельно учить weighted mean/SD и ridge coefficients; lambda1/10/100 выбирать по pooled weighted inner-heldout Brier, ties1e−12→больший lambda. Затем заново выучить preprocessing и final fit на всех ранних primary rows.

Late outcomes и late weights не должны входить в fit signature, нормализацию, tuning, fit/preprocessing, stopping tolerance, выбор модели или выбор domain. Отдельная prediction-функция получает late features и ранний artifact. Временная late validation матрицы допустима; её значения не должны менять ранние model parameters. В имеющемся начальном `forward32/model.py` `fit_forward` не принимает late outcomes/weights; inner tuning и preprocessing используют early arrays, `predict_forward` имеет только fixed-artifact scoring. Это положительная структурная проверка, ещё не полный аудит будущего pipeline.

Необходимы синтетические тесты: изменение late features не меняет `artifact.models`/early training metadata; изменение late outcomes на уровне orchestrator не меняет ранние параметры; scoring отдельно не мутирует artifact; late labels не передаются через общий DataFrame. `late_domain.n` может законно меняться с размером scoring matrix, поэтому тест инвариантности относится к fitted models, а не к всему envelope artifact. Train/test person-key disjointness и masked-PSU disjointness проверять внутри runner, не публикуя ключи.

Калибровка на late — **диагностика**, не recalibration: никакого нового intercept/slope для изменения уже scored probabilities, никакой коррекции различий белков по late распределению и никакого late-outcome выбора penalty.

## 4. Два разных режима неопределённости

Оба режима заранее фиксируют500 attempted replicates и отдельно экспортируют все500 строк со статусами; интервалы не объединяют.

| Режим | Что остаётся фиксированным | Что resample | Интерпретация |
|---|---|---|---|
| fixed-trained | Все4point early fitted artifacts и late predictions | Полный late eligible design, внутри каждой stratum m−1 PSU с replacement, multiplier count·m/(m−1), затем фиксированный primary-mask | Поздняя survey-sampling неопределённость **условно на одном наблюдавшемся раннем fit** |
| joint-refit | Feature forms/domain/tuning grid/early inner partition | Независимый full-early-design draw с повторным preprocessing/tuning/final fit; тот же late draw, что в fixed mode | Неопределённость всей зафиксированной ранней обучающей процедуры плюс late sampling, при наблюдаемых эпохах/измерениях и публичном design |

Одинаковые late replicate weights использовать для всех4arms **и обоих режимов**, а train draws получать отдельным RNGstream/seed. Это уменьшает различия Monte Carlo при сравнении режимов, но не делает их одной статистической выборкой и не разрешает pooling intervals. Hash folds должны оставаться прежними при bootstrap, включая нулевые replicate weights; не перераспределять surviving PSU между фолдами.

При m=3 реально выбрать m−1=2 PSU draws; нельзя применять shortcut для двухPSU ко всем strata. Draw на eligible frame → primary-mask, а не наоборот. Ни при неудаче fit, ни при неидентифицируемой калибровке нельзя перерисовывать bootstrap sample или менять penalty/fold rule. Сохранять failed attempt и определённость каждой пары/метрики. Предпочтительны arm-local statuses; если orchestration заранее делает joint draw all-or-nothing при ошибке одного arm, это отдельно записать и не выдавать за четыре независимо неудавшихся fits.

Percentile bounds — приближённые design-bootstrap интервалы, а не доказанная finite-sample95% coverage. Здесь nonlinear penalized fit, дискретный выбор lambda, ranking/ties AUC, ограниченное число поздних masked units и возможные separation/слабая определённость calibration. Joint bootstrap учитывает повторное обучение, но не выбор feature forms/маркера по предыдущим исследованиям, неизвестную assay calibration, хранение, новые будущие эпохи или межкогортный transport.

## 5. Поздние события, calibration и интерпретация

Новые late n/events в этом review не просматривались. Ранее pooled общая выборка не гарантирует точность отдельной поздней оценки. Перед score export отдельно показать unweighted events/non-events, weighted observed risk, n/effective n, positive-weight events и число strata/PSU. Не считать большое weighted population count большим числом независимых событий.

Primary Brier математически определён при ненулевой сумме весов даже в single-class draw; AUC требует positive-weight event и non-event. Calibration slope/intercept может быть неидентифицируемой при separation, почти постоянном score или неудаче optimizer. Не превращать такие diagnostics в нулевые значения и не подавлять valid/undefined counts. Все40 заданных metric/contrast interval targets в каждом режиме должны иметь собственные valid/undefined. Если bounds построены лишь на определённых draws, это conditional defined-draw bounds, не полный nominal95% CI. Нулевая частота ошибок прогнозного fit не означает500 определённых calibration slopes.

До метрик зафиксировать математические stop/status rules, а не после просмотра событий. Не менять горизонт/популяцию или начинать объединять периоды из-за широкого late CI. Пустой early/late primary domain, неизвестный source/schema, нарушение8y uncensoring или невозможность обязательного early fit — blocker соответствующего empirical path; редкие late события сами по себе ограничивают вывод, а не разрешают post-hoc редизайн.

## 6. Assay/routine drift не устраняется временным split

Сохраняются определения31: cycle1999 serum correction `1.013*raw+0.147`, identity2001/2003, только serum; creatinine-only CKD-EPI2021, UACR из urine reported units, исходный HbA1c без withdrawn regression, неизменный research SBP mean и smoking coding. Форма compact11 и четыре arms не меняются из-за late метрик. Первичные основания записаны в `b2m31/MEASUREMENT_REVIEW.md` и `CREATININE_REVIEW.md`.

Для SSCARD_A известно: baseline serum1999–2004, измерение2018–2020 на SiemensVista, mg/L, flags0/1/2, исправленныйJuly2022релиз. Не установлены release-specific lot/IFCC traceability, распределение исходных циклов по assay batch или sample-level QC; flag0 — detection-range status, не доказанный LOQ/QC-pass. Длительность хранения различается с baseline calendar cycle; без batch/storage metadata нельзя разделить biological epoch drift, storage effect и laboratory drift. Одинаковая платформа/эпоха assays не доказывает межцикловую аналитическую эквивалентность.

Следовательно допустим лишь reported-values forward prognostic analysis с явной оговоркой. Не переносить коэффициентыIII/старых cystatin assays, не рассчитывать Cys-derived eGFR, не называть это трёхбелковой панелью: BTP здесь отсутствует. Возможное ухудшение late performance не доказывает биологическую бесполезность маркера; возможное улучшение не доказывает assay transport или клиническую готовность.

## 7. Что необходимо закрепить и экспортировать

1. Commit/LOCK до сравнительных расчётов: cycles, target/common00,8y endpoint, все4arms/2contrasts, transforms, native weights, fixed3fold hash/tuning, оба bootstrap modes/seeds,500attempts/status rules, запрет recalibration и post-hoc subgroup/partition search.
2. Четыре final early artifacts: ordered features, learned means/scales/constant flags, coefficients/intercepts, lambda grid/inner scores/tie choice, convergence gradients и inner-fold aggregate support; source/plan/code fingerprints. Никаких person vectors/IDs.
3. Source-audit aggregates32: фактические full-release SSCARD detection/pristine/weight marginals, native W4/W2 positive/missing counts и weight sums по циклам и partition; общий/reported domain flow, пересечения причин exclusions, weighted coverage/effective n; отдельный full→eligible masked inventory каждого периода; два mortality readers и8y gate. Эти наблюдённые агрегаты устраняют предел31, где whole-release marginals имелись лишь как runtime guard.
4. Late point metrics/paired contrasts/calibration bins; wide bootstrap aggregates для **обоих отдельных modes**, узкие primary/secondary CSV и независимый stdlib interval/status/hash audit. Полные500 попыток и metric-local failure counts обязательны; не только успешно fitted draws.
5. Прежние frozen models, результаты и научные файлы не менять. Повторные запросы двум denied mortalityPDF не нужны: допустимо перенести явно документированный amendment31 с прежними primary definitions/hashes/two-reader gate, сохранив отдельные статусы `denied_not_retried` без фальшивого SHA/placeholder.

Поддерживающий независимый review весов: `/tmp/forward32_weights_review.md` (подготовлен другим агентом по уже имеющимся primary metadata, без GET/данных/метрик).

**Итоговый gate:** дизайн можно фиксировать и выполнять после перечисленных спецификаций. Нового доступа/аккаунта/соглашения не требуется по известному публичному пути. На данном этапе нет denied/pending/uncertain новой операции: новые сетевые или raw-действия в этом review не предпринимались. Исторические дваPDF имеют прежний исход **denied HTTP403, not retried**, не текущий необработанный запрос.
