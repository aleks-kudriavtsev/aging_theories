# Итерация34: готовность к независимой репликации и условная точность

**Независимая репликация пока не выполнена:** среди18 конкретных releases/visits
не квалифицирован ни один открытый независимый person-level источник для
routine±B2M±Cys и underlying cancer002/heart001. Это ограниченный аудит названных
маршрутов, не утверждение об отсутствии данных во всём мире. Публикационные
таблицы и codebook не использованы как participant validation. Нейтральный33B
не изменён и не переинтерпретирован; split-search, новые маркеры/исходы и новые
модельные fits не запускались.

Реализованы pure source-readiness gate, conditional paired-loss precision tool,
12 явно условных сценариев,4 bounded prediction-atom numerical checks, source/
code manifests и frozen method proposal с нерешёнными условиями активации.
Лабораторный reportv3/презентация не затрагивались; Fine–Gray/joint отложены.

## Что проверено в источниках

Начальная карта —freeze30; теперь BTP не нужен. Для **algorithm replication**
с новым обучением иной документированный клинический concentration assay сам
по себе не blocker. Не перенесены старые CEP/brand/full3 gates. Для неизменных
коэффициентов отдельно нужен source scorer и assay/unit mapping;40OOFfits33B
не один finalfrozenmodel. OlinkNPX/SomaRFU/DNAm не подменяют mg/L measuredpair.

Public001/002/010 coarsening достаточно для двух targets; known-other010 —
competingnoncase, unknowncause — не010. Нужны broadheart,malignantunderlyingcause,
не CHD alone/allCVD/contributingcause/cancerincidence. Source-specific baseline,
same-visit routine11,96month endpoint/censoring,design/weights/family clusters
проверяются раздельно. Разные visits/assayyear не создают новых участников.

|Route|Release/visit|Access action /grant|Blockers|Verdict|
|---|---|---|---:|---|
|NHANESIII_FIRST|1988-1994 first exam with stored-serum releases|not_requested /not_required_for_anonymous_public_release|10|not qualified|
|NHANESIII_REPEAT|1988-1994 repeat examination SSNH3CYS/37a|not_requested /not_required_for_anonymous_public_release|12|not qualified|
|NHANES_SSCARD1999_2004|corrected July2022 SSCARD_A biological baseline1999-2004|not_requested /not_required_for_anonymous_public_release|1|not qualified|
|NHANES2005_2014|2005-2014 public catalogues/2019LMF|not_requested /not_required_for_anonymous_public_release|12|not qualified|
|NHANES2015_2018|2015-2018 public catalogues/2019LMF|not_requested /not_required_for_anonymous_public_release|12|not qualified|
|ARIC_V2|ARIC visit2 approximately1990-1992 /BioLINCC2026a|not_requested /not_granted|15|not qualified|
|ARIC_V4|ARIC visit4 approximately1996-1998 /BioLINCC2026a|not_requested /not_granted|15|not qualified|
|ARIC_V5|ARIC visit5 approximately2011-2013 /BioLINCC2026a|not_requested /not_granted|15|not qualified|
|ARIC_V6|ARIC visit6 2016-2017 /BioLINCC2026a events through2023|not_requested /not_granted|15|not qualified|
|ARIC_V7|ARIC visit7 2018-2019 /BioLINCC2026a events through2023|not_requested /not_granted|15|not qualified|
|CHS1992_1993|CHS1992-1993 protein visit /2026BioLINCC package|not_requested /not_granted|15|not qualified|
|MESA_KIDNEY2012_2014|late kidney ancillary2012-2014 /standard events2020|not_requested /not_granted|15|not qualified|
|FHS_EXAM7|Offspring exam7 measured-protein lead /BioLINCC events2021;Omni visit unverified|not_requested /not_granted|15|not qualified|
|UKB2006_2010|baseline2006-2010 clinicalCys and OlinkB2M metadata|not_granted /not_granted|14|not qualified|
|HRS2016_VBS|2016VBS /retained30access+trackerreview|prohibited /prohibited|15|not qualified|
|AGES_KIDNEY2010|683-person2010-2011 concentration study,distinct fromlargeSomaScan baseline|not_requested /not_granted|15|not qualified|
|CRIC_BASELINE|retained30clinicalkidney-marker publication;currentNIDDK package unverified|not_verified /not_verified|17|not qualified|
|UCI_CKD|public400-rowCKDclassification24-feature schema|not_requested /not_required_for_anonymous_public_release|13|not qualified|

Каждый blocker с verified_false/not_verified, точной причиной и sourceIDs сохранён
в READINESS_RESULTS.json и CANDIDATES.json; ни одно unknown не объявлено absent
вообще. SOURCE_MANIFEST.json содержит32URLs,7review hashes и provenance старого
каталога. Это ledger прочитанных/retained metadata, **не rawsource snapshots**:
originalsourceSHA=null, originalbytes не сохранялись. Два ранее deniedmortalityPDF
и oversized/reCAPTCHA пути не обходились; новых person-level GET не было.

Конкретные пути для пользователя, **не выполненные этой работой**:

- ARIC/CHS/MESA/FHS/AGES: eligible humanPI/учреждение должны проверить текущую
  eligibility/цель, получить exactrelease/proposal+IRB/EC+approval и требуемый
  RMDA/DMDA/DUA. Статусnot_requested/not_granted, неpendingapplication и неtool
  denial. Затем нужны actualsamevisitpair/core и underlyingcause/time/schema.
- ARICV2:UACR unqualified;V4:HbA1c unqualified;V5:B2M не подтверждён LongLab;
  V6/V7:baseline2016–19 слишком поздний для8лет до2023cut. Не подменять core
  glucose/старым визитом или followup другого baseline.
- CHS: consentCVD-only не разрешает автоматически cancer; currentcontinuousage/
  BMI и samevisitHbA1c/UACR неqualified. Legacybins не заменять midpoint.
- MESA-Kidney2012–14: exactlateancillary,alignedcore и более поздний endpoint
  vintage нужны отдельно от standardevents2020. FHS: measuredOffspringexam7
  plasma lead существует, но units/QC/jointrelease/cause остаютсяunqualified;
  Omni visit не объявленexam7. Pedigrees требуют собственногоclustercontract.
- UKB: нет предъявленного approvedproject/processing scope; clinicalCys30720
  существует, B2M225впрочитанномOlink —NPX, неqualifiedclinicalmg/L. Подтверждение
  доступа не исправляет assay mismatch. Текущую возможность новой заявки
  проверяет humanapplicant; прежнееpauseнаблюдение не выдаётся за новыйдопуск.
- HRS: providerAI/LLM prohibition остаётся hardrouteblocker. Новыйlogin/DUA не
  разрешает этот person-levelAIroute; возможен только фактически разрешённый
  human-led маршрут вне запрета сallowedaggregates. Никаких обходов.
- CRIC: currentterms/release неизвестны после прежнегоmetadata403, который здесь
  не retry; humanreview exactrelease/erratum/assays/endpoint нужен преждеэмпирики.
- Открытые NHANESIII/repeat/SSCARD уже использованы; разрешение не превратит их
  вuntouchedcohort. УUCI_CKD нет нужнойprotein/mortalitypanel; большееN этого
  классификационногонабора не создаёт отсутствующихисходов.

Ключевые primaryroutes: [BioLINCC FAQ](https://biolincc.nhlbi.nih.gov/faq/),
[ARIC2026a](https://biolincc.nhlbi.nih.gov/studies/aric/),
[CHS](https://biolincc.nhlbi.nih.gov/studies/chs/),
[MESA](https://biolincc.nhlbi.nih.gov/studies/mesa/),
[FHS](https://biolincc.nhlbi.nih.gov/studies/framoffspring/),
[UKB clinicalCys](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=30720),
[UKB NPX](https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=1839).
Точные supportingdetails/границы —ACCESS_REVIEW/ASSAY_REVIEW иpre-resultaddendum.

## Что именно рассчитывает precision tool

Дляfixed prediction rules D=(Y−p_new)^2−(Y−p_ref)^2 в[-1;1]. При заданных
pairedvariancev,total lossDEFF иevaluableN:
SE≈sqrt(DEFF*v/N), halfwidth≈z*SE,
N_conditional=ceil(z²*DEFF*v/h_target²). Effectsign иincidence сами по себе не
задаютv илиN. AnchoredSE требует N_ref/DEFF_ref; marginalscoreSD безcovariance
не принимаются. ОбщийDEFF уже включает заявленные источникиinflation: не
умножать его ещё раз наKishcorrection. Методprecision разных performance
measures должен задаваться отдельно ([Riley2021](https://pubmed.ncbi.nlm.nih.gov/34031906/));
pairedmeanformula здесь — элементарная заявленнаяvarianceapproximation, не
приписываемая этой статье клинически validatedextension.

**Для primaryalgorithmreplication достаточностьN остаётсяunqualified**: обучение,
подборλ иoverlappingtrainingfolds не описываются varianceindividualD/N. Это
conditionalproxy, не полный training-aware samplesize. AUC/calibrationprecision
изBrierне выводятся. User/stakeholders должны заранее выбрать осмысленныецели;
в34 никакой clinicalthreshold/miminaleffect от собственного имени не задан.
Все значения ниже —illustrativetechnical assumptions, не33B-derivedeffects/SEs,
не новый empiricalrisk и не мощность33B. LOW/HIGH —sensitivityrows, неCI/priors.

|Scenario|AssumedμD|AssumedSDD|DEFFloss|Technicalh|CurrentN_eval|Assumedrawπ|ConditionaltargetN|Recruitment expectation|Events atCURRENT N_eval /predictive range|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
|CENTRAL_TECHNICAL_NULL|0|0.02000|1.5|0.0005|5000|0.030|9220|12294|150 /[127; 174]|
|LOW_RAW_EVENT_PROBABILITY|0|0.02000|1.5|0.0005|5000|0.010|9220|12294|50 /[37; 64]|
|HIGH_RAW_EVENT_PROBABILITY|0|0.02000|1.5|0.0005|5000|0.050|9220|12294|250 /[220; 281]|
|HYPOTHETICAL_BENEFIT|-0.0005|0.02000|1.5|0.0005|5000|0.030|9220|12294|150 /[127; 174]|
|HYPOTHETICAL_HARM|0.0005|0.02000|1.5|0.0005|5000|0.030|9220|12294|150 /[127; 174]|
|LOW_LOSS_SD_DESIGN|0|0.01000|1|0.0005|5000|0.030|1537|2050|150 /[127; 174]|
|HIGH_LOSS_SD_DESIGN|0|0.05000|3|0.0005|5000|0.030|115244|153659|150 /not applicable: dependent events|
|LOW_JOINT_ASSAY_CORE_YIELD|0|0.02000|1.5|0.0005|5000|0.030|9220|23050|150 /[127; 174]|
|LARGER_EVALUABLE_DOMAIN|0|0.02000|1.5|0.0005|20000|0.030|9220|12294|600 /[553; 648]|
|WIDER_ILLUSTRATIVE_HALFWIDTH|0|0.02000|1.5|0.001|5000|0.030|2305|3074|150 /[127; 174]|
|EARLY_CENSORING_METHOD_BLOCK|0|0.02000|1.5|0.0005|5000|0.030|9220|13660|150 /[127; 174]|
|ANCHORED_SE_CONDITIONAL_FIXED_RULE|0|0.03162|1.5|0.0005|5000|0.030|23049|30732|150 /[127; 174]|

Eventranges относятся к явно показанному **CurrentN_eval**, не к ConditionaltargetN. Это5000 дляобычныхстрок и20000 дляlargerdomain. Чтобы оценить counts при другомN, нужно явно задатьэтоN, не читать соседнююtargetколонку какrawtrials.

Например, SD_D=.02,DEFF1.5,h=.0005 дают9220evaluablepairs; q_joint=.75 и
q_known=1 дают12294candidates **по ожиданию**. Это не рекомендуемыйN для
рефит-модели или гарантия достижения выборки. Acrossspecifiedassumptions
conditionalN1537–115244 —диапазон допущений, неconfidenceinterval. Нулевой,
положительный иотрицательныйhypotheticaleffect при том жеv даютодинаковоеN.

q_joint —совместная доступностьcore+обоих assays; произведение marginalrates
безнезависимости недопустимо. q_known —endpointknownnessconditionalonthatpanel;
оно уже учитываетendpointloss. Деленияна(1−c)второйразнет. Даже10%раннего
censoring выводятIPCW_method_extension_required=True, а неготовуюCIFvalidation.
Unknowncause/неизвестныйstatus не превращаются вnoncase приувеличенииN.

Expectedraw events=N*pi сpi **вsampleevaluablepopulation**, неweightednational
incidence. Binomialpredictiveranges условны поindependentBernoulli events, не
incidenceCIs/power. LossDEFF≠eventDEFF: dependent-eventcase не получает exact
binomialrange, и N/DEFF не объявляетсяactualtrials. Иллюстративныйeventfloor100
вinputтолькоcount-taildiagnostic; ни100,нипрежние50events неpowercertificate.

## Numerical validation: bounded synthetic laws, не новые cohorts

Предсказания заданы доY наfixedXatoms; затемY~Bernoulli(r). ExactD moments
включают within-Yvariance иvarianceconditionalmeans. Sharedclusterpairmixture
имеет exactDEFF=1+(m−1)rho. По20000реплик наcase, фиксированныеseeds20261120–23.
Все12предзаданныхmomentchecks (meanD,varianceDbar,meanevents) прошли6MCSEband;
coverageнеиспользовалоськаккритерийуспеха иseed/случаи не выбирались по нему.

|Toy law|N /clusters|ExactDEFF|Known-variance normalCI coverage|EstimatedclusterSE normalCI coverage|MC95%Wilson range for estimatedSE coverage|
|---|---:|---:|---:|---:|---|
|IID_VALID_PREDICTION_PAIRS|1000 /1000|1.00|94.525%|94.470%|[94.145%; 94.778%]|
|EXPLICIT_POSITIVE_CLUSTER_MIXTURE|1000 /200|2.00|94.965%|93.685%|[93.339%; 94.014%]|
|RARE_EVENT_SKEW_DIAGNOSTIC|1000 /1000|1.00|96.225%|86.915%|[86.441%; 87.375%]|
|IDENTICAL_RULE_DEGENERATE|100 /100|1.00|100.000%|100.000%|[99.981%; 100.000%]|

Вrare-eventcase nominal95%normalinterval соценённымSE покрывает истиннуюmean
только86.915% (MCrange86.441–87.375%). Это честный результат заданного toyDGP,
не неисправность, замаскированная сменойseed/threshold. Он показывает предел
normalplanningapproximation, **не оценивает покрытиеbootstrap33B** и не объясняет
его нейтральныйрезультат через post-hocpower. DegenerateD0 даёт100%mathematical
coverage, а не доказательство превосходства95%метода. Фактическиеassays,
реальнаяsurveydependence иrefittedalgorithmprecision этим не валидированы.

## Воспроизводимость исходного расчёта и изолированная публикация

Protocol4df636ea18ca8bb29b05788c37800438566426ec; scientificcode
b054c89b42d6a73f5ed5caf6d8ccd0e0f0b8430b. Локальный productionpackage выполнен
один раз, Python3.12.14/NumPy2.3.5/SciPy1.17.0; пять exactstdoutJSONblocks
сохранены. UnitpushCI40936928789:41tests eachPython3.11/3.13 success.
**Scenariosjob в этомCI skipped**, поскольку он PR-only. Не называть локальный
package запущенным вCI и не утверждать scientificCI длянего наfinaldocshead.
107scientificfiles/103priorhashes сохранены; старые33Bmetrics/артефакты неизменны.
Независимый локальный аудит260/260PASS: пять exactstdoutfiles,12формул/PMFquantiles,4atom/cluster/MCSE/Wilson laws,двухрежимные18sourcegates иsource/codehashes. MC/модели аудитор не пересчитывал и новыеoutcomes не читал.


Эти результаты, `RESULT_PROVENANCE.json` и аудит260/260 сохранены **побайтно
как исторический локальный расчёт**. Их прежний authorization blocker описывает
состояние того запуска. Последующая явная авторизация разрешила отдельный PR
пакета34 в `aging_theories` и его слияние после проверок. Текущий статус и CI
изолированной публикации записываются отдельно в `PUBLICATION_PROVENANCE.json`.
Ветка `aging_biomarkers` в эту авторизацию не входит.

Прямая ветка исследований отличалась отmain на437файлов; этот пакет включает
только итерацию34. Старые модели и данные **не включены вmain**. Для22исторических
metadata/code pins и103прежних Pythonфайлов нужен отдельный checkout точного
коммита27ea33fc6cfd594c309b63ae59e695129511a9c8. Runner читает и проверяет его
хеши, но не импортирует его модели. Расчётные входы, seeds, transforms и gates
34 не изменены; адаптер только явно разрешает пути к историческим зависимостям.
CI повторяет сценарии и сравнивает четыре научных/source JSON с сохранёнными
оригиналами побайтно. Для summary допускаются только документированные версии
окружения и хеши кода изолированной упаковки. Это повторная проверка того же
расчёта, не новая эмпирическая репликация.

Из checkout изолированного пакета:

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0
git fetch origin 27ea33fc6cfd594c309b63ae59e695129511a9c8
git worktree add --detach ../replication34_history 27ea33fc6cfd594c309b63ae59e695129511a9c8
# Исторический checkout должен иметь HEAD точно27ea33fc...; его файлы проверяются SHA256.
export REP34_HISTORY_ROOT="$(cd ../replication34_history && pwd)"
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_replication34*.py' -v
OPENBLAS_NUM_THREADS=1 python -m research.mortality.replication34.analyze --history-root "$REP34_HISTORY_ROOT" --out /tmp/new_replication34_results
python tools/verify_replication34_publication.py --history-root "$REP34_HISTORY_ROOT" --results-dir /tmp/new_replication34_results --out /tmp/new_replication34_publication_audit.json
```

`tools/verify_replication34_aggregates.py` — неизменный аудитор исходного
исследовательского checkout661bcaf94e59a8c6a18bda0cc9771c7515b83989. Его
сохранённый результат260/260 не подменён аудитом публикации; он проверяет
исходную топологию исследовательской ветки, которой вmain нет.

Frozenmethodproposal —REPLICATION_PROTOCOL.md. Этоblockedproposal/template,
не cohort-specificpreregistration: qualifiedrelease,mappings/design ивыбранные
precisiongoals должны бытьзафиксированы ДОновыхoutcomes. Primarymode —algorithm
replication; fixed-sourceкоэффициентныйtransport отдельнотребуетfinal8scorers,
которых в34необучали. Никакихactualexternalpredictions,modelupdates илиclinical
readinessclaims. Только metadata,assumptions,aggregatehypotheticalsimulations.

Следующий разрешённый содержательный шаг — либо фактически разрешённый exact
independentrelease с закрытием указанныхgates, либо отдельно предописанная полная
generative simulation обученияалгоритма дляоценки его traininguncertainty. Не
новыйпереборстарыхразбиений/маркероврадиположительногорезультата.
