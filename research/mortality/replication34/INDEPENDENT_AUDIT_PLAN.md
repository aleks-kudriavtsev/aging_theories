# Итерация34: независимый методический план и план аудита

Подготовлен2026-10-29 **до новых результатов34**. Цель этапа — квалификация независимой репликации и прозрачные сценарии точности. Это не новый подбор на33B и не оценка эффекта по искусственным данным. Новые person-level данные, predictive outcomes/metrics, raw/ZIP GET, refits и Fine–Gray не используются. Labreportv3, презентация, main и прежние модели вне scope.

На момент подготовки каталог `research/mortality/replication34` и immutable протокол ещё не появились. Сверка с конкретными JSON-контрактами, кодом, final scientific commit и CI остаётся pending; этот документ не объявляет их проверенными. Прочитаны retained freeze30 qualification/source reviews. Отдельный агент подготовил `/tmp/rep34_precision_audit_expectations.md` и8 независимых арифметических fixtures без эмпирических входов.

## 1. Разделить readiness и расчётную точность

Qualification отвечает: существует ли конкретный разрешённый независимый release с интерпретируемыми совместными Cys+B2M, routine core и нужными исходами? Precision отвечает: какова ожидаемая ширина интервала при **заданных гипотетических** моментах paired loss и design/availability assumptions? Даже идеальный сценарий precision не открывает доступ, не доказывает пригодность assay и не является результатом внешней репликации.

Текущая панель — **Cys+B2M**, BTP не требуется. Нельзя механически использовать прежний отрицательный full3 вывод freeze30 как готовый ответ для двух белков. Freeze30 остаётся неизменным источником проверенной карты и конкретных access/assay/endpoint ограничений; gate34 должен заново применить именно текущие требования.

Нужно назвать вид будущей репликации: воспроизведение заранее фиксированного marker contrast с новым обучением на иной когорте либо применение замороженного scorer. Это разные требования к assay bridge и параметрам. Существование40OOF fits33B само по себе не даёт единственного внешнего frozen scorer. До выбора и фиксации операции нельзя объявлять deployment/transport readiness.

## 2. Fail-closed qualification gate

У каждой обязательной категории — отдельные статус `confirmed / blocked / unknown`, primary URL/точный release, дата и объяснение. Unknown не преобразуется в pass. Overall empirical-ready false, пока хотя бы один обязательный gate не подтверждён; допускаются отдельные readiness-to-protocol/metadata статусы с чёткой семантикой.

| Gate | Что необходимо доказать |
|---|---|
| Access и разрешённая среда | Person-level release действительно доступен в текущем разрешённом режиме; open article/codebook не равен open microdata. Нет новых аккаунтов, подписания terms, заявок, restricted data или контактов. AI/person-level prohibition сохраняется; пользовательский boolean не отменяет политику держателя. |
| Independent provenance | Иной recruitment/release и отсутствие прежнего использования в программе; учёт NHANESIII иContinuousNHANES, повторных визитов, surplus re-assays иподвыборок. Новое имя/assay year не означает новых людей. Не пытаться идентифицировать/связывать людей ради проверки overlap. Distinct sampling frame не доказывает буквально нулевой возможный person overlap; предел metadata evidence раскрывается. |
| Joint baseline panel | Cys иB2M у **тех же людей** ипри совместимой исходной точке, с известными matrix/units/assay/release/QC/LOD правилами. Marginal n двух assays иразных визитов не дают joint n. NPX/RFU не называть mg/L; не придумывать bridge/IFCC traceability. |
| Routine mapping | Возраст/пол, smoking,SBP,BMI,HbA1c, serum creatinine с source-specific calibration иcreatinine-only eGFR, urine albumin/creatinine с правильными единицами. Не подставлять glucose/CHD proxies/прошлый визит безотдельного заранее названного core ипротокола. |
| Time иcause endpoints | Exam/baseline origin,96месяцев, vital status, administrative censoring иunderlying cause с обоснованным mapping cancer/heart. Не считать allCVD эквивалентом heart001 или incidence cancer эквивалентом cancer death002; MCOD не underlying cause. Unknown early cause ираннее цензурирование требуют stop либо отдельного заранее заданного estimator/protocol. |
| Domain/design | Конкретная популяция, общий joint panel/core/outcome domain, weights/strata/PSU или иной sampling design; наблюдаемый event support. Нельзя считать50событий универсальной достаточностью. Возрастной/болезненный domain shift раскрывается и требует соответствующего estimand. |

Готовность не объявляется по календарному запасу follow-up или по одному максимальному death date. Cohort gate не должен выгружать records для проверки доступности. Если все доступные пути blocked/unknown, empirical branch имеет конкретный terminal blocker; stage34 всё равно может закончить qualification infrastructure и **явно гипотетические** precision calculations.

## 3. Парный Brier contrast: математические ограничения

Определить знак один раз:

`D = (Y − p_expanded)^2 − (Y − p_reference)^2`.

ПриY∈{0,1} иp∈[0,1]: **D∈[−1,1]**. Поэтому допустимые assumed meanμ иSDσ удовлетворяют `|μ|≤1`, `σ≥0`, `σ²≤1−μ²`. Эти bounds — mathematical feasibility, не эмпирическая оценка. При|μ|=1 допустима толькоσ=0. Не допускать impossible moments, NaN, negative n/SD или fractions вне[0,1]. Boundary cases либо обрабатываются точно и помечаются degenerate, либо явно отклоняются контрактом; нельзя молча получитьnan/inf или «идеальную» реальную precision.

Pairing существенно: нужна SD **разности** losses, а не две независимые SD или их сумма. Incidence сама по себе не определяет эту SD, correlation прогнозов, calibration илиразмер эффекта. Нельзя брать старую33B оценку effect/SD/SE как скрытый/default вход. Все scenario assumptions имеют явный источник `illustrative / external design assumption / explicitly supplied`; ни один не выдаётся за новое observed effect.

## 4. Аналитическая precision и двойной учёт design

Для выбранного marginal95% normal approximation:

`z = NormalDist().inv_cdf(.975) ≈1.959963984540054`

`SE ≈ σ * sqrt(DEFF_total / n_evaluable)`

`halfwidth ≈ z * SE`

`target_n_evaluable ≈ ceil(z² * σ² * DEFF_total / h²)`.

Это **assumption-based approximate precision**, не guaranteed recruitment, power или рекомендация клинического sample size. Effectμ не входит в fixed absolute halfwidth формулу; он влияет лишь на bounded feasibility илина заранее заданный relative-h target. Не вычислять post-hocpower, chance ofsignificance илиN дляподтверждения старого33B знака.

DEFF должен иметь определённый scope. Еслиinput ужеdesign-adjusted SE илиeffective n, нельзя повторно применятьDEFF/Kish/clustering. При SE-input нужны его reference n иdesign scope; нельзя масштабировать его на инойn безфиксированной scaling assumption. Unequal-weight иcluster corrections нельзя перемножать безоснования, еслиtotalDEFF уже включает оба. Stratification может даватьDEFF<1: simulator обязан либо иметь подходящий DGP, либо честно не поддерживать этотслучай; нельзя незаметно округлятьDEFF кцелому или использоватьn/DEFF какточное число наблюдений.

Normal CI может выходить за допустимый parameter domain; не клиппировать его толькоради вида. Если заранее задана intersection с[−1,1], показывать исходный интервал/метод и не объявлять это доказательством номинального покрытия. Small event count/число независимых кластеров иnon-normal paired losses остаются ограничениями, даже приформально маломhalfwidth.

## 5. Availability, censoring ичисло событий

`q_joint_assay_core` — именно joint fraction, не произведение неизвестно зависимых assay/core marginals. Если дополнительно используется `q_endpoint_known`, он должен быть **условным наjoint panel/core domain**. Тогда `E[n_evaluable] = n_recruited*q_joint*q_known` — scenario accounting, не гарантия realised n и не поправка selection bias. Еслиfraction0, цельнедостижима вданномсценарии; это не повод делитьна0 иливключать synthetic records.

`E[events]=E[n_evaluable]*p_case`, `E[noncases]=E[n_evaluable]*(1−p_case)` толькопри явно заявленной sample-event probability данного retained domain. Survey-weighted population incidence не всегда является sample-event fraction; оба понятия не подменяются. Если сценарий одновременно описывает две взаимоисключающие причины в однойпопуляции, их incidence sum не может превышать1; отдельные сценарии не объявляются joint model.

Known competing early death — известный noncase, **не потерянный follow-up**. Censoring/unknown-cause fraction может описывать ожидаемый отбор в endpoint-known subset, но простоеумножение на(1−c) не доказывает unbiased8yCIF и не заменяет IPCW/survival estimator. Если реальная новаякогорта имеетраннеецензурирование, нынешний complete96m binary protocol нельзя просто применить послеудалениялюдей: gate останавливает dependent empirical analysis либо требуется новый estimator/protocol. Это различие должно быть явным ивschema, ивотчёте.

## 6. Валидация распределения/симуляции

Симуляция — synthetic planning diagnostic, не marker validation. Требуются фиксированные seed,number ofreplicates,явные assumptions,allattempted outputs иno seedsearch. Не использовать clippedGaussian как будто это заданное bounded Brier-difference distribution: clipping меняет mean/SD иможет искусственно улучшить precision.

Допустимый пример математического DGP: `D=2*Beta(a,b)−1`, где m=(μ+1)/2,v=σ²/4,k=m(1−m)/v−1,a=mk,b=(1−m)k. Дляinterior case0<v<m(1−m) обаshape>0;σ0 deterministic,maximumvariance — отдельный endpoints-only Bernoulli boundary. Есливыбран другой bounded DGP, егоsupport иmoments проверяются аналитически.

DEFF>1 требуетобоснованной dependence/cluster construction либо явно обозначенной analyticinflation, не скрытого количества «effective iid records». Actual cluster sizes/ICC→DEFF relation иrounding должны быть проверяемы. Monte Carlo mean/variance сравниваются спредписанными analytic moments сзаранее заданной MC tolerance; normal-CI coverage проверяется тольковусловияхданного DGP иCIprocedure. Coverage, еслирассчитана, имеет собственный MonteCarloSE/interval и не является post-hocpower. Согласие sample meanSD сформулой не валидирует calibration/assays/clinicalmodel илиналичие independent cohort.

## 7. Независимые sanity fixtures

- Увеличениеn в4раза припрочихравных вдвоеуменьшаетSE; DEFF в2раза увеличиваетSE вsqrt(2);halfwidth target вдвое меньше требует примерно4×n доокругления.
- Приσ=.02,DEFF1,h=.001: unrounded target≈1536.5835283,ceil1537. ПриDEFF2:ceil3074. Эти числа — только algebra fixtures,непредложениеразмераисследования.
- Приn_recruited10000,q_joint.5,q_known.8,p_case.02: expected evaluable4000,expected events80,noncases3920. Это не гарантия80реальныхсобытий и не certificate adequatepower.
- Разныеμ приодинаковыхσ,DEFF,h даютодинаковый absolute-precision targetn, еслиmoments feasible.
- μ=.6,σ=.8 — variance boundary;μ=.6,σ>.8 impossible. μ=±1,σ>0 impossible. Fraction0/1,incidence0/1,invalidnegative/bool inputs,degenerateσ0 иn1 имеютявную семантику.
- Математически identical input scenarios воспроизводятся byte/structure-wise; перестановка cohort labels илиcase labels не может менять precision arithmetic, еслиnumeric assumptions одинаковы.

## 8. Финальный аудит после протокола

После появленияimmutable replication34 protocol/LOCK читаютсяони, код исценарии. Затем независимый stdlib auditor проверяет:

1. Pin сохранённых103scientific modules/33B результатов иfrozenmodels; exactscenario inputs доoutput; новыеgates отдельно отfreeze30. BTP нигде неостаётся обязательным current-panel criterion.
2. Каждое cohort decision восстанавливается изmandatory statuses,access/assay/core/endpoint/overlap evidence; impossible/unknown input не превращаетсявready. Fresh primarysource citations идатаотделеныотretained30 metadata; log/raw retrievallimits названыточно.
3. Всеanalytic SE/halfwidth/inverse-n,availability/event accounting,rounding иbounds пересчитываютсяотдельнойmath; nonidentified assumptions/SE scope раскрыты.
4. Simulation DGPparameters/support/moments,seed/replicatecount/DEFF construction иMC diagnostics проверяютсянезависимо; simulation output нельзяпредставитькакреальнуюcause marker benefit.
5. Exactcompleted scientificCI plaintextblocks→repo bytes/SHA; synthetictests иCI привязаныкscientificcodehead,а неlaterdocs-onlyhead. Никакихraw/ZIP/personrecords/refits audit; толькоmetadata/aggregates идопустимые synthetic calculations.

Итоговыедоставки: `/tmp/rep34_audit.json/md` иstandalone verifier послесогласованиясхемы; донихэтотметодическийплан не заявляет новыхрезультатов илиготовойрепликации.
