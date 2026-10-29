# Cause33: source/endpoint review для отдельных восьмилетних исходов

Дата: 2026-10-29. Scope: только документация и исходный код; без person-level
файлов, ZIP, новых cause counts, model fits, сравнений или изменений репозитория.
Основание кода: aging_theories `fe4d20672dc761c8a6533b3dcaf5a6a0340251d8`
(сохранённая итерация32), неизменённые source/measurement gates31.
Совместная модель причин и Fine–Gray отложены: здесь они не реализованы и не
запускались. Прежний остановленный joint-аудит не является новым протоколом.

## Вывод

В SSCARD_A baseline1999–2004 технически доступны девять именованных групп
первоначальной причины смерти, из которых можно заранее определённым правилом
выбрать **не более пяти**. Код010 — неоднородный остаток, а пропуск причины —
отдельная неопределённость: они не являются дополнительными именованными
заболеваниями. Фактический top5 и event support в нужном общем домене пока НЕ
установлены этим review; показатели иной выборки/горизонта их не заменяют.

При полном наблюдении восьмилетнего статуса и известной причине смерти отдельный
бинарный исход `I(dead=1, t_exam<=96, cause=k)` имеет среднее, равное
восьмилетней cumulative incidence соответствующей **опубликованной** группы.
Известная смерть от другой причины до96месяцев — известный noncase этого
исхода, не цензурирование. Это допускает отдельно названную бинарную модель
риска на фиксированный горизонт. Её коэффициент не является cause-specific
hazard ratio; вероятности отдельных моделей нельзя складывать или объявлять
совместным распределением причин. Старые all-cause baselines31/32 сохраняются.

## Первичные источники и границы проверки

1. [NCHS2019 public LMF dictionary,28.04.2022](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf).
   Группы/ограничение по baseline подтверждены текстом официальной страницы,
   возвращённым поисковым индексом CDC в текущем read-only metadata review.
   Сам PDF не открывался/не загружался заново; checksum/визуальная проверка
   не заявляются. Ранее отклонённый запрос PDF остаётся `denied_not_retried`.
2. [Текущая официальная NCHS страница linked mortality](https://www.cdc.gov/nchs/data-linkage/mortality-public.htm).
   Раздел2019public подтверждает публичную версию и изменение follow-up/cause
   для части записей. Vital status при этом не изменялся, но это по-прежнему
   результат linkage с обозначениями assumed alive/deceased. Текущий2022
   linkage — restricted; его наличие не добавляет открытые исходы и не
   разрешает смешивать vintages.
3. [Официальный SSCARD_A codebook, исправлениеJuly2022](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.htm).
   Baseline1999–2004; assay2018–2020; SSB2M/SSCYST serum mg/L. Веса
   WTSSCB4Y/WTSSCB2Y относятся к этому исправленному subsample. Флаги
   detection не являются индивидуальным QC/LOQ certificate. Это тот же
   ранее использованный биологический материал, не новая внешняя когорта.
4. [Официальный Weighting Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx).
   Таблица1999–2004 и примечание о subsample weights поддерживают формулуW6
   ниже; для1999–2002 требуется предоставленный четырёхлетний вес.
5. [Официальный Variance Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx).
   Design задаётся masked `(SDMVSTRA,SDMVPSU)`; доменные оценки и
   неопределённость не равнозначны simple-random-sample анализу. Репликационная
   исследовательская схема не становится официальным NCHS variance estimator
   только вследствие группировки по PSU.

Два ранее отклонённых mortality PDF не запрашивались повторно и не заменялись
обходным URL. Description PDF отдельно не переоткрывался. Доступ/аккаунт/
соглашение/контакт с держателями не запрашивались.

## Точные опубликованные группы

`UCOD_LEADING` описывает первоначальную underlying cause, а не contributing
cause. Кодбуковый тип — трёхсимвольный код; текущий проверенный reader хранит
его численно как1.0…10.0. Для новых target IDs безопасно канонизировать
валидированное целое к `001`…`010`, сохраняя исходное значение.

| UCOD | Именованная опубликованная группа | Ограничение интерпретации |
|---|---|---|
|001|Болезни сердца|Не вся cardiovascular mortality и не отдельно ИБС|
|002|Злокачественные новообразования|Нет локализации опухоли|
|003|Хронические болезни нижних дыхательных путей|Не отдельноJ44/ХОБЛ|
|004|Несчастные случаи|Не все внешние причины, не отдельный суицид|
|005|Цереброваскулярные болезни|Нет ишемического/геморрагического подтипа|
|006|Болезнь Альцгеймера|Не вся деменция|
|007|Сахарный диабет|Underlying cause; не MCOD-флагDIABETES|
|008|Грипп и пневмония|Совместная группа, не все инфекции|
|009|Нефрит, нефротический синдром и нефроз|Широкая группа, не отдельноХБП/N18|
|010|Все остальные причины|Остаток; не отдельная defined disease дляtop5|

Точные ICD/recode отношения уже заданы в `endpoint_contract11.py`; здесь не
вводятся новые подтипы или замены широких кодов узкими диагнозами. В частности,
старые пять heads12 —001,002,003,005 и остаток004/006–010/unknown — не означают,
что в новом вопросе можно ранжировать только эти пять или считать остаток
именованным заболеванием.

### Применимость по baseline, не по календарному году смерти

| Baseline NHANES | Public2019 UCOD labels | Пригодность для этой задачи |
|---|---|---|
|1999–2014|001…010|Различимы девять именованных групп и остаток|
|1999–2004 SSCARD_A|001…010|Нужная marker/cohort комбинация; фактические counts надо квалифицировать|
|2015–2016|Только001/002/010|003/005/009 и прочие скрыты в остатке, не нулевые события|
|2017–2018|Только001/002/010|То же coarsening;2019follow-up также не даёт всем живым8лет|

Ограничение относится к году исходного обследования. Смерть в2017–2019 у
участника baseline1999–2004 не теряет остальные коды только из-за года смерти.
Нет основания восстанавливать скрытый subtype по national shares, biomarker
уровню, MCOD-флагам или совпадению public/restricted записей.

## Что уже сохраняет код31/32

В `validate_temporal08.independent_mortality_reader` используются поля
`elig`15, `dead`16, `cause`17–19, `t_exam`46–48 (позиции1-based). Пропуск
`''`/`.` становитсяNaN; `001` становится1.0. Все mortality поля, включая
`cause`, соединяются1:1 с DEMO в `b2m31.cohort.load` и остаются в frame,
возвращаемом `forward32.cohort.load`. Новый download для появления cause не
нужен. Это проверка кода, не новое чтение записей.

При этом прежний all-cause gate проверяет eligibility/status/time, но не
является новым cause-specific semantic gate. Дляcause33 необходим отдельный
строгий валидатор кодов и status/cause согласованности; наличие числовой
колонки само по себе не доказывает готовность конкретного исхода.

Public2019 source содержит disclosure substitutions follow-up/cause для части
записей. Следовательно, оценка касается опубликованного coarsened/perturbed
endpoint; точная индивидуальная причина/дата и ошибка этой замены не
восстанавливаются. Её неопределённость не покрывается обычным bootstrap.

## Fail-closed qualification до marker comparisons

1. Зафиксировать протокол, точную популяцию, ranking rule, horizon96months,
   source/metadata/code SHA и список001…009 ДО новых comparative metrics.
   Не менять старые source/model/results. Не повторять два deniedPDF запроса.
2. Сохранить общий feature-domain31: screening age40–75, examination/linkage,
   11routine terms, positive reported Cys/B2M и оба flags0. Для всех arms и
   всех cause comparisons — одинаковые строки и folds. Любое изменение домена
   требует явного нового estimand/protocol, а не тихого cause-specific drop.
3. Требовать `ELIGSTAT=1`, `MORTSTAT` строго0/1 и finite `PERMTH_EXM>=0`.
   У assumed-alive должен быть follow-up>=96. Смерть с0месяцев остаётся
   событием; нельзя сдвигать её на0.5месяца ради удобства binary model.
4. На этих baseline допустимы только целые001…010 или задокументированный
   пропуск. Код0/11/дробный/бесконечный блокирует branch. Assigned cause у
   alive — конфликт, а не дополнительный event. Новые данные не исправлять
   по предполагаемому распределению исходов.
5. До modelling экспортировать только агрегаты: all-cause events<=96,
   counts/weighted masses001…010 и unknown, по каждому baselinecycle и
   выбранному train/evaluation part. Проверить суммы ordinary/weighted
   causes+unknown против all-cause. Никаких participant IDs/prediction exports.
6. **Неизвестная причина смерти<=96 не является known noncase.** Для точного
   сравнения underlying-cause CIF такой случай блокирует зависимый расчёт,
   пока не зафиксирован честный иной estimand/handling. Не объявлятьNaN
   кодом010, живым или цензурированным. Отдельно можно описать bounds нагрузки
   при отнесении unknown к каждой группе, без вымышленных личных labels.
   Unknown cause смерти после96 не мешает известному noncase на8лет.
7. Для pooled1999–2004 burden использовать именноW6 изSSCАRD:
   `(2/3)*WTSSCB4Y` для1999/2001 и `(1/3)*WTSSCB2Y` для2003. Не брать MEC
   weight, не нормализовать каждый цикл/причину отдельно. **Frame32 имеет
   `weight`=nativeW4early/W2late:** его нельзя выдать за pooledW6; брать
   сохранённый исходныйW6 или заново явно строить формулу из source fields.
   Если сохраняется forward split32, trainingW4/evaluationW2 остаются разными
   native-period estimands, и их burden не надо складывать как единыйW6.
8. Ранжирование в явно объявленном общем шестилетнем домене:
   `M_k=sum(W6*I(dead<=96,cause=k))`, `risk_k=M_k/sum(W6)`; denominator общий
   для всех причин, включая survivors/known competitors. Raw event countD_k
   показывается рядом, tie rule фиксируется заранее. Код010 и unknown не
   входят в девять кандидатов. Показать все девять ranks и причины blocked,
   а не только удачные белковые результаты. До фактического агрегатного gate
   никакие пять причин или их ранги здесь не утверждаются.
9. Если top5 выбирается по pooled outcome burden, это outcome-informed
   выбор исследуемых outcomes в уже использованной когорте. Зафиксировать
   ranking stage до marker metrics и назвать последующий анализ exploratory;
   обычный CI не учитывает выбор ранга/причины. Для заявления о prospective
   forward evaluation список причин должен выбираться только в early training
   с nativeW4; поздняя full-domain таблица burden тогда только описательная.
10. Draft33A уже фиксирует более строгий minimum event-count gate до counts:
   >=50positive-weight target events и>=50known noncases в общем домене,
   event-bearing>=20maskedPSUs; каждый outer evaluation>=5/5,
   outer training>=20/20, inner training>=10/10, inner evaluation>=3/3
   events/noncases. Это computational/reporting thresholds, не power/EPV/
   calibration certificates. Показать все failed cells и концентрацию событий.
   Нехватка support означает blocked cause, не biologically useless marker и
   не повод выбирать favourable split. Порог не ослаблять после counts.
11. Конечные predictions оцениваются без перекалибровки. AUC требует оба
   positive-weight класса; calibration может быть неопределена из-за малого
   support/flatness/separation/сходимости. Во всех resampling попытках хранить
   paired failures и metric-local valid/undefined; defined-only percentile
   bounds не выдавать за полные nominal95%CI. All-cause accuracy не заменяет
   cause-specific support.

## Допустимые отдельные marker contrasts после qualification

Для каждого заранее выбранного допустимогоk первичный marker contrast:
`routine+B2M minus routine`, одинаковые11routine terms иlog2reportedB2M.
Secondary: `routine+Cys+B2M minus routine+Cys`, только в том же заранее
определённом общем домене. Нельзя выбирать основной contrast после просмотра
результатов, ранжировать белки как победителей или объявлять отсутствие
secondary CI сигнала общей бесполезностьюB2M.

Для полного восьмилетнего наблюдения target-кодировка такова:

| Опубликованный статус к96месяцам | BinaryY_k |
|---|---|
|Смерть с известнымcause=k к96включительно|1|
|Смерть с другой известнойcause к96|0, competing noncase|
|Живой сfollow-up>=96 либо смерть после96|0|
|Смерть<=96 сunknowncause|Неопределённо; dependent gate blocks|
|Живой сfollow-up<96|Неопределённо; dependent gate blocks|

Это следует непосредственно из определения `F_k(96|X)=P(T<=96,J=k|X)`;
обучаемая binary expectation равна этому риску при известных labels.
Независимые по исходу regressions могут дать сумму рисков>1, поэтому только
отдельная оценка Brier/AUC/observed/predicted/O-E/calibration для каждогоk.
Не использовать их сумму как all-cause score, survival или residual. Прямой
восьмилетний binary risk — отдельная модельная семья, не прежний cause-specific
engine, не subdistribution-hazard fit и не новая внешняя валидацияfrozen30.

## Статус

Metadata/code feasibility — положительная для отдельногоcause gate в раннем
SSCARD_A. Фактические named-cause counts/ranking/unknown support — ещё НЕ
проверены в этом review; это следующий разрешённый этап после фиксации протокола.
Отдельные cause models, comparative metrics, causal interpretations и clinical
readiness здесь не получены и не заявлены. Результаты32/all-cause сохранены.

## Read-only review draft33A перед первым protocol commit

Прочитаны `research/mortality/cause33/GATE_PLAN.md` и
`CAUSE_DEFINITIONS.json`. Критических endpoint/source противоречий не найдено.
Все десять записей definitions точно совпадают с существующим
`endpoint_contract11.PUBLIC_CAUSES` поid/label/ucod113/icd10; ranking разрешён
только001…009. В плане корректны unknown-early блокировка, baseline-зависимое
coarsening, pooledW6 вместо native32weight, top-five-before-support без замены
удобной шестой причиной, shared domain и отдельный будущий comparison protocol.

Проверены только кодовые определения fold salts: `b2m31:outer` и
`b2m31:inner:{fold}` соответствуют исходникам31 и его прежнему плану. Численные
fold assignments, source records и новые counts не читались/не вычислялись.
Gate-only scope честно отделён от marker benefit inference. Тексты primary
limitations выше следует сохранить в итоговом gate report, особенно
perturbed causes/times и отсутствие нового доступа к двум deniedPDF.
