# Итерация34: доступ и endpoint readiness для ARIC, CHS и MESA

Проверено **2026-10-29**, только первичные официальные страницы и metadata.
Прочитаны существующие freeze30 reviews; broad catalogue заново не составлялся.
Assay-часть согласована с `/tmp/rep34_assay_review.md`. Записи участников,
ZIP и новые исходы не запрашивались; нет аккаунтов, принятых terms, заявок,
соглашений, контактов с держателями или model fits. Репозиторий не изменён.

## Решение для текущего контракта

В трёх проверенных близких маршрутах **не найден уже разрешённый анонимный
person-level package**, закрывающий одновременно routine11, клинические
количественные Cys+B2M, исходную точку и два полных восьмилетних underlying-cause
исхода cancer002/heart001. Это terminal blocker немедленной новой empirical
репликации в текущем режиме. Это не глобальное доказательство отсутствия
подходящей когорты и не отрицательный результат полезности белков.

**BTP не требуется.** Старый трёхбелковый freeze30 gate нельзя использовать
как готовый отказ текущему двухбелковому исследованию. Но удаление BTP не
создаёт согласованный визит, более длинный cause follow-up или доступ.

Основной режим34 — **репликация алгоритма с новым обучением**: заранее те же
11routine terms, четыре arms, два B2M-контраста и training-only nested tuning
в действительно иной когорте. Это отличается от scoring неизменными33B
коэффициентами. У33B сохранены40point OOF fit records, а не один заранее
определённый final-fit scorer. Нельзя задним числом выбрать/усреднить их и
назвать это замороженной моделью для внешней проверки.

Для нового refit иной изготовитель концентрационного assay сам по себе не
является отказом: новый training учит свои means/SD/coefficients. Нужны
интерпретируемые измерения, материал/единицы, метод, calibration/QC и policy
reported values с честными ограничениями. Количественный межметодный bridge
особенно существенен для frozen-coefficient scoring. RFU/NPX, DNAm-прокси и
другие относительные сигналы — отдельный вопрос, не клинические mg/L по умолчанию.
NHANES-specific CEP/1999 serum corrections, веса и masked PSU labels нельзя
переносить в новую когорту без её собственных source definitions.

## Доступ: почему «Open BioLINCC Study» не означает открытые записи

[Текущий BioLINCC FAQ](https://biolincc.nhlbi.nih.gov/faq/) различает открытое
чтение study pages/dictionaries и запрос данных. Для последнего нужны
institutional-email Login.gov, допустимый PI, protocol, IRB/EC determination,
одобрение NHLBI и RMDA до release. Для некоторых проектов достаточно exempt
review, но это не отменяет документацию и approval. FAQ также отсылает к
действующим ограничениям NIH CADR по местонахождению/учреждению; их применимость
не устанавливается по имени или языку пользователя. Одобрение цели не является
автоматическим разрешением передавать microdata в любую AI/cloud среду.

Текущий статус всех трёх person-level маршрутов: **not_requested /
not_granted**, не pending application и не tool-level denial. Metadata pages
прочитаны успешно. Ни «free/no cost», ни публичная статья, ни открытый codebook
не являются разрешением на anonymous participant download.

## Фокусная матрица актуальных releases

| Маршрут | Что подтверждено текущей официальной страницей | Что нужно разрешить до расчёта |
|---|---|---|
| **ARIC — HLB00020026a** | [BioLINCC](https://biolincc.nhlbi.nih.gov/studies/aric/), updated03Aug2026: visits1–10, annual follow-up2–35, mortality/heart/stroke events through2023, ancillary studies; Request Data требует login. Commercial use запрещён. Исходный recruitment45–64лет, четыре communities, n15792 — это не joint routine/protein n. | Конкретный same-visit core и pair; underlying cause/death/censor поля именно выданной версии; доступ/цель/среда. Более поздние B2M visits не получают автоматически8лет по наличию долгой общей истории ARIC. |
| **CHS — HLB00040025a** | [BioLINCC](https://biolincc.nhlbi.nih.gov/studies/chs/), updated23Sep2026: follow-up throughyear35, но **events through30Jun2015**; recruitment65+, n5888. Consent tiers различают CVD-related-only и unrestricted research; commercial restrictions также есть. | Cancer scope нельзя предположить для CVD-only consent. Exact routine/pair/continuous-age/cause release не квалифицированы; два разных follow-up cutoffs нельзя слить в один «35лет событий». |
| **MESA — HLB00640825a** | [BioLINCC](https://biolincc.nhlbi.nih.gov/studies/mesa/), updated15Dec2025: exams1–6, events through2020. Среди20 перечисленных ancillaries есть Cys#057; kidney/B2M joint package прямо не установлен. Request Data требует login; commercial investigator consent tiered. | Отличать original baseline от MESA-Kidney2012–14, подтвердить inclusion именно нужного ancillary, aligned clinical visit, exact cause coding и censor cutoff. Наличие Cys ancillary не доказывает B2M или joint11-domain. |

### ARIC: сильный measured-pair lead, но пока не полный контракт

Согласованный assay-review по [официальному LongLab](https://www5.cscc.unc.edu/aric9/sites/default/files/public/visitdocuments/v11/ARIC%20V1_V11_LongLab_261210_2_np%20Derived%20Variable%20Dictionary.pdf)
установил B2M mg/L и Cys наV2/V4/V6/V7. V5 B2M в этих строках не перечислен:
это **unverified**, не доказательство отсутствия любого V5 ancillary. V4
same-visit HbA1c, V2 UACR не квалифицированы; нельзя заменить HbA1c glucose или
взять мочу с другого визита под прежним routine11.

[Официальный V7 Manual17](https://www5.cscc.unc.edu/aric9/sites/default/files/public/visitdocuments/v7/ARIC%20MOP17%20190206.pdf)
задаёт2018–19; [Manual30,p81](https://aric.cscc.unc.edu/aric9/sites/default/files/public/visitdocuments/v9/ARIC_Manual_30.pdf)
указывает началаV6=2016/V7=2018. Поэтому current endpoint cutoff2023 не даёт
живым этих visits полный96месячный исход. Увеличивать его воображаемыми более
новыми outcomes нельзя.

Найдена конкретная endpoint lead: [официальный INCBY24 derived dictionary,
§3.6UCOD](https://aric.cscc.unc.edu/aric9/sites/default/files/public/surveillancedata/Data%20Dictionary%20INCBY24%20Cohort%20Surv%20Incident%20Derived.pdf)
описывает underlying cause из death-certificate field, с NDI fallback.
Это подтверждает study-level конструкцию, **не подтверждает**, что exactUCOD,
date и censor fields включены в конкретный BioLINCC2026a release. Полный11.2MB
main/18MB ancillary dictionary повторно не запрашивался: прежний metadata
reader size limit не обходился.

**Конкретное будущее действие пользователя:** допустимый PI выбирает точный
визит/release и запрашивает absolute Cys/B2M, все same-visit routine inputs,
underlying ICD cause, death/last-contact/admin-censor dates и соответствующий
consent/design documentation. ДляV6/7 отдельно нужен достаточно поздний
фактически существующий cause/censor vintage. Пока такого approved package
нет, empirical path остановлен; новый account/request здесь не создавался.

### CHS: consent, текущая кодировка и baseline alignment — отдельные gates

[Direct CHS distribution policy](https://chs-nhlbi.org/CHS_DistribPolicy)
требует association со study; collaborator получает sponsor, paper proposal,
approval, institutional IRB и DMDA. Все фактические users должны быть указаны.
Это отдельный маршрут от BioLINCC, а не anonymous bypass.

Согласованный assay-review не установил exact same-visit clinical mg/L
Cys+B2M release/QC; присутствие Cys плюс SomaScan target этого не доказывает.
Для protein baseline1992–93 same-visit HbA1c/UACR не квалифицированы.
[Официальный MOP index](https://biolincc.nhlbi.nih.gov/media/studies/chs/Manual_of_Operations.htm)
перечисляет Year9SpotUrine: это предостережение против подстановки позднего
UACR в ранний baseline, не доказательство отсутствия иной ранней ancillary мочи.

Есть дополнительная source-specific неопределённость:
[legacy public value labels,pp1–6](https://biolincc.nhlbi.nih.gov/media/studies/chs/hidden/labels/valuelbl.pdf),
до сих пор связанные с current study через
[Other Documents](https://biolincc.nhlbi.nih.gov/media/studies/chs/Other_Documents.htm),
описывают age в двухлетних категориях и collapsed BMI/anthropometric tails.
Это документация старого public release, **не утверждение**, что нынешний
approved package одинаков. Точные continuous age/BMI нужно проверить;
midpoint возраста или «восстановление» хвостов не являются source values.

Exact underlying heart/cancer mapping нынешнего package здесь не установлен.
Обнаруженное в indexed main dictionary поле Cause of Death само по себе не
показывает его полный coding/coverage. FatalCHD, CVDdeath и cancer history
не заменяют два требуемых underlying mortality endpoints.

**Конкретное будущее действие пользователя:** approved PI получает совместимый
с cancer consent tier, exact assay/core visit и death/censor/cause package;
проверяет текущую privacy coding возраста/BMI и sampling weights. При выборе
direct CHS route нужны sponsor/proposal/IRB/DMDA, при выборе BioLINCC — его
approval/RMDA. Ни один путь не инициирован.

### MESA: baseline2000 и kidney2012 — разные исходные точки

[Официальный overview](https://mesa-nhlbi.org/about/overview) указывает original
baseline2000–2002, recruitment45–84 без knownCVD, наблюдение mortality и
получение death certificates из state vital records/NDI. Это потенциальная
независимая cohort lineage; n6814 не является нужным paired-biomarker n.

В coordinated assay-review [первичные kidney methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5572311/)
подтверждают реальные концентрации Cys и B2M, но для поздней MESA-Kidney
подвыборки2012–14. Same-visit routine11 и inclusion нужных файлов в current
standard package не квалифицированы. До2020 от поздних kidney visits нет
полных8лет у всех живых; follow-up нельзя начинать раньше даты крови.

Indexed official [MESA2025a dictionary](https://biolincc.nhlbi.nih.gov/media/studies/mesa/data_dictionary/MESA_2025a.pdf)
показывает `dthtype`, `dthtext`, `dthtt` и время до смерти/lastfollow-up.
Полный code mapping cancer002/broadheart001 в данном review **не подтверждён**.
История сбора death certificates — сильнее aggregate paper, но слабее точной
released endpoint specification. Полный15.3MB main dictionary здесь не
загружался: использованы явно ограниченные indexed metadata, не полный review
выдаваемого package. Alternate raw-file канал не использовался.

**Конкретное будущее действие пользователя:** допустимый PI запрашивает именно
согласованный kidney/clinical/death release, а не лишь общий MESA package;
подтверждает endpoint cutoff достаточный от фактического kidney baseline,
ICD/cause definitions и applicable consent. Альтернатива — документированный
другой same-visit pair, если он существует; здесь его наличие не установлено.

## Общие scientific gates для разрешённой будущей выдачи

1. **Независимость.** Отдельный recruitment/release не означает автоматически
   нулевое пересечение лиц с NHANES или другими cohort studies. Здесь cross-study
   IDs не сверялись и re-identification не предпринималась. Нужны provenance и
   история использования; новый визит тех же людей не считается второй когортой.
   В изученном roadmap ARIC/CHS/MESA были metadata candidates, не завершёнными
   empirical tests программы; это не exhaustive claim обо всех частных доступах.
2. **Популяция и time zero.** Age определяется на выбранном specimen/clinical
   baseline. ARIC recruitment age, CHS65+ и MESA45–84 не равны автоматически
   target40–75. Нужен объявленный overlap domain и число retained participants;
   extrapolation за75 или перенос более поздних routine measurements не скрываются.
3. **Два endpoint.** Нужна underlying cause, с заранее проверенным ICD9/10
   crosswalk к malignantC00–C97 и heartI00–I09/I11/I13/I20–I51. CHD-only уже,
   а CVD соstroke шире heart001; cancer incidence/history не cancer death.
   Достаточно двух сохранённых target groups и корректного known-other;
   **все девять именованных причин из gate33A здесь не обязательны**.
   В частности, publicNHANES2015–2018 сохраняет001/002/010: это coarsening
   само по себе не блокирует cancer/heartCIF. Код010 — known competing other,
   а не unknowncause. Для этих поздних циклов другие отдельные барьеры —
   полный8y public2019follow-up, joint clinical pair и история использования.
   Unknown early cause не становится known noncase. Полный96месячный binary
   outcome требует no early-censored survivors; иначе нужен отдельный protocol
   изменения estimator, а не удаление коротко наблюдавшихся.
4. **Routine и assay.** Same-visit age/sex/smoking/SBP/BMI/HbA1c/serum creatinine/
   urine albumin+creatinine и измеренные clinical Cys/B2M с понятными единицами.
   Новые методы допустимы для нового refit при достаточной интерпретируемости;
   отсутствие per-personQC certificate не заменяется вымышленным pass flag.
   Неподтверждённые поля остаются unresolved, а не объявляются отсутствующими
   во всей когорте. Ни glucose/history proxy, ни поздний UACR, ни NPX/RFU не
   подставляются под прежний контракт.
5. **Дизайн и веса.** Нужна cohort-specific policy: biospecimen/subsample
   selection, weighting, clusters/семьи/centres и целевой received-domain
   estimand. NHANESW6 и87maskedPSUs не имеют силы в иной когорте. Unit weights
   тоже должны быть обоснованы заранее; нельзя выдумать национальную
   репрезентативность, pseudo-PSUs или family independence.
6. **Precision до comparison.** После разрешённой выдачи сначала source/domain/
   event-support gate: target events/known noncases, event-bearing groups,
   weight concentration, missing/censor/cause ambiguity. Общий размер cohort
   и опубликованные mortality totals не заменяют joint routine+pair+cause n.
   Затем immutable comparison protocol/commit до новых predictive metrics,
   одинаковые people/folds во всех arms и полный учёт uncertainty/undefinedness.

## Что выполнено и точный blocker

Выполнена focused current-release/access/endpoint metadata проверка трёх
маршрутов, сверка с прежним roadmap и согласование нового двухбелкового
algorithm-replication mode. Public pages, metadata dictionaries и первичные
публикации не классифицированы как person-level anonymous files.

Блокируемое действие — получить и обработать новый независимый joint
routine+Cys+B2M+heart/cancer microdata package для8y репликации. Цель —
внешняя cohort replication алгоритма и precision, а не split search на33B.
Причина — нет предъявленного approved release/operator/environment и остаются
перечисленные same-visit/endpoint/horizon проверки. Статус **not_requested /
not_granted**; новых person-level действий с исходомdenied/pending/uncertain
не было. Необходимое следующее человеческое действие описано для каждого
маршрута; агент его не выполнял и не принимал terms. Если такое разрешение
недоступно, эта empirical ветка терминально заблокирована в текущем режиме.

33B остаётся нейтральным exploratory результатом; metadata readiness не
превращает его в доказанный выигрыш. Main, старые модели/результаты,
labreportv3/presentation и отложенные joint/FG не затрагивались.
