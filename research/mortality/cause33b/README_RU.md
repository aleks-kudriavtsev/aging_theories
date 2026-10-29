# Итерация33B: отдельные восьмилетние CIF для cancer002 и heart001

**В этой exploratory проверке убедительная добавочная прогностическая ценность
B2M не установлена ни для одного из двух заданных контрастов и исходов.**
Point Brier снижается, но все восемь marginal95% интервалов ΔBrier/ΔAUC
пересекают0. Это не доказательство нулевого эффекта/эквивалентности или общей
бесполезности маркера. Различие эффектов между причинами не тестировалось:
численные размеры/знаки и интервалы не заменяют interaction/contrast.

## Общий домен и корректный исход

Те же6044 участника NHANES1999–2004/SSCARD_A,age40–75,MEC baseline,flag00
complete-case domain; eligible6546,all-cause8y586. Cancer002176events/75eventPSU;
heart001148/67. Horizon96months inclusive. Weighted observed incidence
cancer2.383497%,heart1.882176%. Corrected six-year weights2/3WTSSCB4Y early,
1/3WTSSCB2Y2003–04;43maskedstrata/87PSU включая третьюPSU stratum1.
Сохраняется received complete-case/measurable estimand, не безусловная national
population оценка и не рейтинг России/мира. Один denominator для всех8fits.

Yk=I(dead1,t_exam<=96,UCOD=k). AUC cases — именно эти deaths; controls —
**все известные noncases**, включая ранние конкурирующие смерти, смерти после96
и survivors сfollowup>=96. Controls5868/5896; ранних competing deaths410/438.
Известный competitor не цензурируется и не называется alive. Unknown common
ранняя причина/ранняя survivor censoring блокируют все сравнения; здесь unknown0.

Это отдельные binary fixed-horizon CIF regressions, **не cause-specific
hazards**, Cox или Fine–Gray. Их probabilities не суммируются вall-cause риск,
выживание или joint-state vector. Исторический пятипричинный piecewise
exponential ridge/CIF engine и all-cause baselines не заменены.

## Все четыре заранее заданных arms

|Исход|Модель|Brier|AUC|Средний прогноз|Calibration slope|
|---|---|---:|---:|---:|---:|
|Cancer002|Routine|0.0228147840|0.72764931|2.405172%|0.975367|
|Cancer002|Routine+B2M|0.0227466185|0.73896316|2.411819%|0.957134|
|Cancer002|Routine+Cys|0.0227578234|0.73558274|2.407605%|0.958293|
|Cancer002|Routine+Cys+B2M|0.0227539537|0.73668094|2.404745%|0.941038|
|Heart001|Routine|0.0177215633|0.77735860|1.875821%|0.958803|
|Heart001|Routine+B2M|0.0175107622|0.78670161|1.849977%|0.932186|
|Heart001|Routine+Cys|0.0174634618|0.78622701|1.847039%|0.945969|
|Heart001|Routine+Cys+B2M|0.0174543737|0.78941827|1.864322%|0.968332|

Core11/units/creatinine corrections иlog2reported serum Cys/B2Mmg/L строго31,
без новой стандартизации/калибровки/imputation/clipping или выбора лучшего
протеина. Detection0 не индивидуальный QC/LOQ pass. Нет exclusive marker→cause
назначения. Не сравнивать absolute Brier между разными-prevalence causes как
рейтинг качества моделей.

## Парные контрасты и500 полных refits

|Исход|Контраст|PointΔBrier|Marginal95%CI|PointΔAUC|Marginal95%CI|
|---|---|---:|---|---:|---|
|Cancer002|B2M сверх routine|-0.0000681655|[-0.0002370227; 0.0001668272]|0.01131385|[-0.00997081; 0.03926355]|
|Cancer002|B2M сверх routine+Cys|-0.0000038697|[-0.0001482213; 0.0002113158]|0.00109820|[-0.01678131; 0.01609796]|
|Heart001|B2M сверх routine|-0.0002108011|[-0.0005336795; 0.0001253487]|0.00934301|[-0.00942779; 0.04192884]|
|Heart001|B2M сверх routine+Cys|-0.0000090880|[-0.0003101009; 0.0001754368]|0.00319126|[-0.00764398; 0.02526044]|

Первичный вопрос для каждой причины — B2M сверх routine поweighted8yBrier;
вторичный — B2M сверх routine+Cys. AUC/calibration — diagnostics. Два исхода,
два контраста и несколько показателей: интервалы **маргинальные exploratory**,
не simultaneous/familywise подтверждение, p-values/positive winner selection нет.
Исходы выбраны поfull-cohort33Acounts/burden; selection/programme history не
включены в CI. Все циклы/маркеры уже использованы; это не независимая external
validation, untouchedholdout или перенос frozen30.

500 общих eligible-frame within-stratum whole-PSU draws, seed20261105. Одни
и те же rows/folds/drawweights для2causes×4arms; zero-weight rows сохранены.
В каждой попытке заново5outer×3inner preprocessing,lambda1/10/100 иouterfits;
weighted mean/SD толькоactualtraining, interceptunpenalized,tie1e-12→largerλ.
Логика и PSUhashes31 неизменны. Нет redraw, fixed-prediction CI, recalibration,
вариантов horizon/subgroups/assayflags или splitsearch. В наличии500wide/
2000pairedCSV, sharedweight SHA/moments и40pointouterfit records.

Все Brier/AUC и paired контрасты500/500defined; predictivefit/metricerrors0.
Совпадение preprocessing подтверждено для всех доступных matched cause/arm/
fold fits; здесь всеpoint иbootstrap predictivefits доступны. На разных
causes могут отличаться lambda/coefficients при одинаковыхX/w/transforms.

## Калибровка без изменения прогнозов

|Исход|Модель|Offset intercept|Full intercept|Slope|Full intercept/slope valid|Offset valid|
|---|---|---:|---:|---:|---:|---:|
|Cancer002|Routine|-0.009474|-0.092715|0.975367|498/500|500/500|
|Cancer002|Routine+B2M|-0.012429|-0.155332|0.957134|497/500|500/500|
|Cancer002|Routine+Cys|-0.010572|-0.150084|0.958293|495/500|500/500|
|Cancer002|Routine+Cys+B2M|-0.009343|-0.206010|0.941038|499/500|500/500|
|Heart001|Routine|0.003660|-0.137435|0.958803|496/500|500/500|
|Heart001|Routine+B2M|0.018928|-0.208934|0.932186|497/500|500/500|
|Heart001|Routine+Cys|0.020608|-0.162132|0.945969|499/500|500/500|
|Heart001|Routine+Cys+B2M|0.010436|-0.096462|0.968332|498/500|500/500|

Все point diagnostics вычислены; quantilebins40строк. Вbootstrap fullcalibration
intercept/slope неопределены1–5/500 наarm из-за несходимости диагностического
полногоfit, offsets500/500. Ошибки сохранены какoffset_only_full_fit_nonconverged,
не predictive failures; соответствующие defined-only bounds условны и не
полные nominal95%CI. Улучшение калибровки не установлено; diagnostic fits
не использованы для исправления OOF scores. Clipping1e-12 — только diagnosticlogit,
не scoredprobabilities. Утверждений clinical readiness нет.

## Протокол, источники и независимый аудит

Protocol **d772b874198a9f7739225bbce5d1b6b6dfcd0c08** до новых прогнозов.
Scientific head **ab09ca0b605b44544f49d21d1161c61609ef92d1**, actualPRmerge
**819b97537bdbf0fec55f50e4aa4a882fd3da67f5**. Scientificrun40928005361:
empirical123638594625 и30tests eachPython3.11/3.13 success; workbench40928005687
success. Автоматический старый33A regression40928005320 такжеsuccess, не новый
model experiment. Последний docs/results-only commit сskipCI не считается новым
scientificCI. Все103scientific hashes равныtestedtree;99старых файлов сохранены.

50official sources/25historical data hashes, unchangedreader/measurement/weight
policy. Вся scientificprojection33A воспроизведена:3974derivedfloats,0differences;
counts/structures/folds/support exact. Preplannedroundofftolerance1e-12/1e-9 не
потребовалась. Ранние deniedmortalityPDF403 не retry; новыхchecksumclaimsдляних
нет. Sourcecause/timeperturbation,assumedlinkage,reportedassay/calibration/storage,
complete-case selection и историческая repeated-use uncertainty сохраняются.
Причины003/004/005 послеfailedgate не fitted; пороги/горизонт/folds не менялись.

Независимый аудит **265/265PASS** сверил6exactstdoutfiles/fileSHA, все80интервалов
из500wide/2000pairedrows,40outerrecords/tuning/preprocessingidentity,control/case
accounting, sources/103codehashes и3974floatprojection. Аудитор отдельно получил
completedGitHubplaintextlogs; rawparticipants не перечитывались, индивидуальные
Brier/AUC и gradients поraw не пересчитывались, refits не повторялись. Это аудит
агрегатов/provenance, не независимая когортная валидация.

ВGitHubchunkedplaintextlog дваU+FEFF стоят передtimestamp; удалены только
transportprefixes, неJSON/CSVcontent. Все5embeddedoutputSHA и6independentexports
побайтно совпали. ZIP не скачивался; artifactAPI digest/expiry записаны отдельно.

```bash
python -m pip install numpy==2.3.5 pandas==2.2.3 scipy==1.17.0 scikit-learn==1.8.0
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_cause33b*.py' -v
python -m research.mortality.cause33b.retrieve --out /tmp/new_cause33b_sources
OPENBLAS_NUM_THREADS=1 python -m research.mortality.cause33b.analyze --source /tmp/new_cause33b_sources --out /tmp/new_cause33b_results
python tools/verify_cause33b_aggregates.py --repo . --log /tmp/cause33B_success_job.log --science-commit ab09ca0b605b44544f49d21d1161c61609ef92d1 --run-id 40928005361 --job-id 123638594625 --tests-per-python 30 --ci-jobs /tmp/cause33B_independent_ci_jobs.json --ci-log-dir /tmp --out /tmp/new_cause33b_audit
```

Дляаудита нужныcompletedjob metadata иplaintextlogs сименами
cause33B_independent_job_JOBID.log. Публикуются только агрегаты, pointparameters,
checksums и код; participantIDs/labels/individualpredictions/weightvectors/rawfiles
не exported. Прежние модели/main/UI сохранены; labreportv3 иpresentation не
затрагивались этойработой. Никаких merge/approvals/deploy/accounts/agreements/
restricteddata/holdercontacts/scheduledtasks.

Следующий содержательный шаг — отдельно зафиксированная проверка чувствительности
обоих тех же контрастов к ограниченному набору PSU-разбиений с публикацией всех
вариантов, без выбора лучшего и без продвижения rarecauses. Это проверит
условность текущих интервалов поединственномуfixedpartition, а не создаст
independentvalidation. Совместная интеграция/Fine–Gray остаются отложенными.
