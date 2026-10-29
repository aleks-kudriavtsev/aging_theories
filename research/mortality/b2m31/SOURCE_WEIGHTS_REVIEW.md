# Iteration31 — квалификация SSCARD_A и весов, metadata only

Дата: 2026-10-29. Записи участников не запрашивались/не читались; изменений репозитория нет. Это подготовка source contract до протокола и сравнительных результатов. Машинное предложение: `/tmp/b2m31_SOURCE_SPEC_PROPOSAL.json`.

## Объединённый шестилетний вес подтверждён

В [официальном NHANES Weighting Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx), Example3 / таблица «Combining three survey cycles», дано объединение 1999–2004 через предоставленный четырёхлетний вес для первых двух циклов и двухлетний для третьего. Примечание непосредственно распространяет эту структуру на subsample weights. Применение к переменным текущего SSCARD_A даёт:

```python
if SDDSRVYR in (1, 2):
    weight_6y = (4 / 6) * WTSSCB4Y
elif SDDSRVYR == 3:
    weight_6y = (2 / 6) * WTSSCB2Y
else:
    raise ValueError("Unexpected cycle")
```

Это обоснованное применение общего официального правила к документированным subsample weights; отдельной готовой шестилетней переменной SSCARD_A нет. `WTSSCB2Y/3` для всех трёх циклов не соответствует правилу: первые два цикла опирались на разные Census population bases. Нельзя брать одноимённые веса из SSBNP_A вместо SSCARD_A без отдельного подтверждения идентичности релиза и домена. После построения W6 не нормализовать каждый цикл отдельно перед общей оценкой: это изменило бы относительные массы циклов. Общая training-only нормализация весов для определённой penalized loss — отдельный выбор протокола, не замена исходного survey weight в оценочных метриках.

Tutorial также требует проверить сопоставимость переменных/популяций и временной тренд при объединении; шестилетний вес сам не решает ни assay bridging, ни missingness. Существующие веса корректируют определённую подвыборку, но не гарантируют отсутствие смещения после дополнительных complete-case/LOD/linked-eligibility ограничений.

## Источник белков и исправление 2022

[SSCARD_A, текущий официальный codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.htm): baseline serum 1999–2004, анализ 2018–2020; eligibility 1+ лет с согласием/сохранённым образцом. Июль2022: удалено 5–<10% ошибочно включённых записей, оставшиеся значения не менялись, веса пересозданы. Нужен именно этот релиз.

Оба белка — Siemens Dimension Vista1500, mg/L. Пределы названы LOD, не LOQ: Cys 0.23–8.00, B2M 0.72–23.0. `*L=0/1/2` обозначает measurable/below/above; ниже LLOD — LLOD/√2. Правило заполнения выше ULOD отдельно не описано. Не превращать опубликованный диапазон в физиологический фильтр или флаг0 в доказательство индивидуального QC.

Вектор кодбука `[flag0, flag1, flag2, missing]`: B2M `[21824,15,24,27]`; Cys `[21824,22,2,42]`. Всего21890; это маргинальные частоты, не совместное взрослое n. `SSPRIS`:20216 pristine,1674 non-pristine. Специальное предупреждение о непригодной non-pristine аналитике относится к glycated albumin; оно не доказывает запрет Cys/B2M. Первичный reported-values/LOD/domain policy и ограниченную pristine sensitivity следует зафиксировать заранее. BTP здесь нет.

Машинный spec хранит эти опубликованные частоты только как проверки полного исходного релиза после разрешённого чтения; это не рассчитанные сейчас персональные агрегаты. Не требовать от будущего общего клинического домена совпадения с маргинальными частотами всего релиза.

## Дизайн выборки — из DEMO

[Variance Estimation Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx) устанавливает public masked variables `SDMVSTRA` и `SDMVPSU` из demographic files. Уникальный кластер задаётся **парой** `(SDMVSTRA, SDMVPSU)`; один SDMVPSU повторяется между strata. Рекомендованный R `nest=TRUE` именно об этом. В SSCARD_A собственных PSU/strata нет; соединять с нужным DEMO по уникальному SEQN и проверять 1:1 / отсутствие unmatched assay keys.

| Cycle | Официальный DEMO codebook | SDDSRVYR | SDMVSTRA, опубликованный диапазон | SDMVPSU, опубликованный диапазон | Полное n DEMO |
|---|---|---:|---:|---:|---:|
|1999–2000|[DEMO](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/DEMO.htm)|1|1–13|**1–3**|9965|
|2001–2002|[DEMO_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/DEMO_B.htm)|2|14–28|1–2|11039|
|2003–2004|[DEMO_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DEMO_C.htm)|3|29–43|1–2|10122|

Не hardcode ровно2 PSU в stratum: в1999 допустим3. Табличные диапазоны сами по себе не определяют фактическое число непустых пар; его надо проверить после получения файлов и повторно в общем домене. Цикл восстанавливается из DEMO/SDDSRVYR, не из даты assay. `RIDSTATR=2` обозначает interview+MEC, а `RIDAGEYR` — возраст на screening, не точный возраст exam; при40–79 он не topcoded. Для иной возрастной популяции нужно пересмотреть определение.

Whole-PSU CV/refit grouping должен сохранять пару во всех вложенных фолдах и при bootstrap-copy. Опубликованный masked PSU не следует называть идентификатором истинной географической PSU. Bootstrap не становится официальным survey-design variance estimator только из-за cluster grouping. Для Taylor domain-оценок tutorial требует сохранять полный survey design и subpopulation indicator вместо преждевременного удаления всех вне-domain строк.

## Mortality и дальнейшие gate

[NCHS public-use2019 LMF description](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf), [dictionary](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf): связь NHANES по SEQN, `ELIGSTAT`, `MORTSTAT`, `PERMTH_EXM`, cutoff31.12.2019. Время/причина некоторых записей perturbированы; vital status не заменён. Общий eligibility/time/8y outcome gate определяется в протоколе и проверяется по данным; календарный запас follow-up не является проверкой фактической полноты исхода. Restricted2022 не используется.

Как установлено в предыдущем read-only аудите, SSCARD_A уже импортирован в `fe414dc1957d8c6041dce48e77cf58ee65b30d28`; Cys оценён в `nhanes04/M4_cystatinC`, B2M загружен, но в предикторной панели benchmark отсутствует. В31 возможен заранее заданный новый B2M contrast на повторной когорте, не новая независимая панель/валидация.

## Direct XPT URL: честный статус

Публичный релиз и download-link presence подтверждены [официальным общим каталогом](https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory), строкаSSCАRD_A/XPT2.2MB, и самостоятельным codebook. Но веб-рендер показывает номер ссылки, не raw href; сохранённого исходного HTML локально не найдено. Существующий benchmark читает локальный filename и не содержит retrieval URL. `source_manifest.json` текущего дерева и setup-копии не содержат SSCARD_A; старый `SOURCE_DATA_SHA256.txt` содержит только hash/filename.

Поэтому spec оставляет `data_url=null`, `direct_link_resolved=false`. Подстановка расширения по шаблону не представлена как проверенный href. Перехода по XPT-ссылке, HEAD/GET к XPT, новых raw-file поисков и повторения прежде отклонённого Python metadata route не было. Прежний исход Python route — **denied**, proxy403; текущий исход raw-запроса — **не предпринимался**, а не pending/uncertain. После фиксации протокола root может разрешённым способом извлечь официальный href из HTML/API перед загрузкой; при сохранении блока источник останется metadata-qualified, retrieval-blocked.

## Проверки, необходимые после протокола

1. Зафиксировать metadata/data URL, фактический hash/bytes обоих; не приписывать старому hash статус новой проверки.
2. Проверить полные размеры/схемы/частоты, уникальные непустые SEQN, cycle membership, допустимые коды flags/PSU; reject неожиданные данные, не пропускать строки.
3. Проверить W4 положителен для SDDSRVYR1/2 и W2 для3; отсутствие W4 у2003–04 не считать маркерным пропуском. Вычислить W6 ровно по формуле и сохранить исходные веса.
4. Отдельно вывести flow, домен/события/пары MVU, веса и missing/LOD/pristine по циклам. Сравниваемые модели используют одних людей и одинаковые фолды; ни один исход не используется для выбора domain.
5. Зафиксировать и выполнить предопределённые reporting/sensitivity policy. Этот документ не выбирает победившую модель и не содержит сравнительных результатов.
