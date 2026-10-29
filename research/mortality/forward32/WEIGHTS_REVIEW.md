# Iteration32: веса и источники для фиксированного forward split

Подготовлено 2026-10-29 **до новых сравнительных метрик**. Использованы только код, закреплённая документация и ранее опубликованные агрегаты итерации31. Записи участников не читались, сетевых запросов к данным, нового обучения и вычисления числа событий по temporal split не было. Это проверка выполнимости и правил, не результат итерации32.

## Правило весов

Для обучения на объединённых **1999–2002** требуется исходный пересмотренный **WTSSCB4Y из SSCARD_A**. Для оценки на отдельном **2003–2004** требуется исходный **WTSSCB2Y из того же релиза**. Официальный [Weighting Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx) требует предоставленный четырёхлетний вес при объединении первых двух циклов из-за перехода Census reference population; правило выбора subsample weights распространяется и на SSCARD. [SSCARD_A codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.htm) определяет эти поля и пересоздание весов в July2022.

Здесь не строится объединённая шестилетняя population estimate. Веса двух частей нельзя описывать как W6: train представляет свой четырёхлетний период, test — свой двухлетний. Не применять WTSSCB2Y/2 к первым двум циклам, не заменять веса на DEMO/MEC или одноимённые поля SSBNP_A. Не использовать label «шестилетняя популяция» для отдельного test score.

```python
if cycle in (1999, 2001):
    temporal_partition = "train_1999_2002"
    analysis_weight = WTSSCB4Y
elif cycle == 2003:
    temporal_partition = "test_2003_2004"
    analysis_weight = WTSSCB2Y
else:
    raise ValueError("Unexpected baseline cycle")
```

Положительность и конечность проверяются только для **нужного этой части** исходного веса. WTSSCB4Y у 2003–2004 ожидаемо отсутствует; одновременный `dropna` по обоим столбцам ошибочно удалит test. При работе с уже возвращённым frame из31 выбранный исходный вес должен быть назначен явно: поле `weight` там содержит W6. `weight6y` можно сохранить для provenance, но нельзя молча оставить его под новым названием исходного четырёх-/двухлетнего веса.

В пределах каждой части W6 отличается от требуемого исходного веса только общим положительным множителем: 2/3 в train, 1/3 в test. Поэтому нормированные weighted averages и training-normalized loss математически инвариантны к такому общему множителю. Это объяснение не заменяет корректные имена/суммы весов и не является новой метрикой. Если penalized loss **не** нормирует веса внутри обучения, масштаб веса меняет относительную силу penalty. Training-only normalization и penalty selection нужно закрепить в протоколе32; test weights не используются для preprocessing/tuning/calibration update.

## Исходные поля сохранены в load31

Статически проверены `b2m31/SOURCE_SPEC.json`, `cohort.py`, `features.py`:

- `SSCARD.required_fields` содержит WTSSCB4Y и WTSSCB2Y.
- `load` соединяет весь выбранный SSCARD frame с DEMO по SEQN и добавляет `weight`/`weight6y`; исходные WTSSCB-поля не удаляются.
- `engineer` начинает с `frame.copy()`, добавляет предикторы и возвращает frame без отбора только feature columns.
- Последний `return d, common, reported, audit, manifest` сохраняет оба исходных веса в d. Сам `load` в этом аудите не запускался.

До новой оценки надо явно проверить наличие исходных столбцов, cycle/SDDSRVYR agreement, нужный положительный вес и ожидаемую пропущенность W4 третьего цикла. При одинаковом domain policy положительное масштабирование не меняет положительность и не должно создавать новую подвыборку; любое отличие counts от ранее закреплённого общего домена требует остановки/объяснения до сравнительных результатов, а не изменения mask по исходам.

## Обязательный дизайн и домен

Masked design берётся из cycle-specific DEMO; ключ кластера — **(SDMVSTRA, SDMVPSU)**. [Variance Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx) объясняет эти public masked variables и nesting. Ни PSU alone, ни SEQN не заменяют ключ кластера.

| Часть | Baseline cycles / SDDSRVYR | Полный публичный masked design из31 | Ранее известный common n |
|---|---|---|---:|
| Train |1999/2001; 1/2|strata1–28, 57 PSU; stratum1 имеет [1,2,3], остальные — по2|4020 = 1889 + 2131|
| Test |2003; 3|strata29–43, 30 PSU, по2 в stratum|2024|

Это **ранее известные** counts/inventory из31, а не новое вычисление событий или оценки модели. Совокупно43strata/87PSU; разбиение непересекающееся и полностью определяется baseline cycle. Число событий и новые temporal scores должны рассчитываться только после protocol commit. Здесь они не выводятся и не восстанавливаются из старых fold exports.

Обязательные gate32:

1. Сохранить тот же target31: screening age40–75, записанный sex1/2, MEC, eligible2019public mortality, корректное exam time и общий finite11-term routine core; оба положительных serum proteins с flags **0/0** для всех arms.
2. Сравнить full DEMO, eligible и common inventory отдельно по частям как множества пар, а не только по числу PSU. Полный/eligible inventory из31 уже совпал; новый temporal adapter должен подтвердить отсутствие ошибочного удаления/переноса пар.
3. Train strata строго1–28, test строго29–43; пересечение пар пусто, объединение равно исходному inventory. Train содержит единственную m3 stratum1; не заставлять все strata иметь ровно2PSU. Singleton или утрата ожидаемой пары — stop, а не post-hoc изменение возраста/периода.
4. Для inner tuning и любых повторных обучений использовать только train PSU и обучение всех преобразований внутри соответствующего train split. Test не участвует в выборе формы, penalty, threshold, split или донастройке калибровки.
5. После protocol commit отдельно вывести flow, missing/LOD/pristine и суммы исходных/analysis weights по циклам/частям. Это позволит проверить применение весов лучше, чем один суммарный denominator31. Вывести фактические события и endpoint gate тогда же, сохраняя96месяцев и отсутствие подмены missing status на alive.

## Assay, storage и ограничения интерпретации

В SSCARD биологический baseline — **1999–2004**, а лабораторные анализы сохранённого serum выполнены в **2018–2020**. Публикация/assay year не определяет temporal partition. Cys/B2M представлены как mg/L на Dimension Vista1500 с analyte-specific flags; это reported measurements. `flag0` означает нахождение в detection range, а не индивидуальный QC-pass, отсутствие storage effect или подтверждение современной standardization.

В полном публичном списке SSCARD variables нет assay batch ID, plate/run ID, индивидуальной даты анализа или storage duration. Поэтому нельзя провести/объявить batch-adjusted forward validation, выводить batch из baseline year или выдумывать переменную assay year для каждого человека. Общее окно2018–2020 не доказывает случайное распределение assay batches между baseline cycles. Долгий срок хранения и возможное смешение cycle, storage и assay ordering остаются неотделимыми по этим release metadata.

SSPRIS отражает pristine/non-pristine specimen; это не batch ID и не точная длительность хранения. Предупреждение codebook о проблемной non-pristine аналитике относится к glycated albumin; оно не доказывает непригодность всех non-pristine Cys/B2M. Primary32 должен сохранять заранее согласованный primary00 domain31; любую storage/pristine sensitivity надо назвать и зафиксировать заранее, а не отбирать по test score.

SSCARD_A и участники уже использованы в31; Cys fitted ещё в прежнем M4, B2M loaded ранее. 2003–2004 использовался и в старом NHANES benchmark, и в31. Поэтому такое разделение — новая **exploratory проверка временной переносимости на повторно используемых данных**, не untouched holdout и не независимая внешняя валидация. Выбор B2M после предыдущих результатов сохраняет selection-history limitation даже при строго forward-only новом fit. Ассайные различия с III и отсутствие BTP не позволяют объявить этот результат валидацией frozen30.

## Итог готовности

Веса, доступность исходных полей и объявленный public masked design **позволяют заранее зафиксировать forward32**. На этой стадии нет нового эмпирического результата. Машинные правила сохранены в `/tmp/forward32_weights_rules.json`; последующие gate должны быть выполнены после protocol commit, без новых source substitutions и без ослабления reader validation.
