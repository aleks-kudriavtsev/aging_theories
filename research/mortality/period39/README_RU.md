# Итерация39: российский возрастно-половой периодный ориентир

Реализован и выполнен строгий adapter и полная таблица для точного возраста40–80,
половых серийF/M и горизонтов5/10лет: **164age/sex/horizon строки**, в каждой
два заранее заданных актуарных расчёта. Это **provisional local-derived period
reference**, а не новая индивидуальная biomarker-модель, disease-free baseline,
текущий календарный прогноз или клиническая российская валидация.

[ГрафикSVG](PERIOD_REFERENCE.svg), [PNG](PERIOD_REFERENCE.png),
[полная таблицаCSV](REFERENCE_TABLE.csv), [машиночитаемые результаты](RESULTS_SUMMARY.json).

## Результат при явно заданных допущениях

Ниже — условные all-cause вероятности synthetic cohort, **если на каждом будущем
возрасте сохраняется schedule2023** из зафиксированного локального производного
CSV. Условие входа — жив на точный день рожденияa. Это не наблюдавшиеся исходы
российской когорты2023→2028/2033 и не вероятность конкретного пациента.
Primary — constant-force approximation внутри каждого годового возраста:

    Q(a,h) = 1 − exp(−sum[m(a), …, m(a+h−1)]).

| Точный стартовый возраст |F,5лет,%|M,5лет,%|F,10лет,%|M,10лет,%|
|---|---:|---:|---:|---:|
|40|1,21|4,11|2,85|9,25|
|50|2,30|6,95|5,47|15,82|
|60|4,80|12,64|11,80|27,42|
|70|11,79|22,80|28,42|46,08|
|80|30,14|40,32|59,85|70,35|

F/M — исходные population sex categories; общеполовая серия не вычисляется
без соответствующих численностей/весов. Для80/10лет используются80…89,
конечный день рождения90. Нецелый возраст не интерполируется. Линии графика
соединяют целочисленные точки лишь для показа.

## Методы и чувствительность

SecondaryUDD сa_x=.5года задаёт однолетнее q=m/(1+.5m), затемQ=1−Π(1−q).
Он не выбран по проверочным исходам; никаких fits/model selection не было.
CF и UDD — **разные** assumptions. Комментарий старого extractor об
«implicitax=.5» для экспоненты не является математическим тождеством; старый
код сохранён, правильная трактовка документирована вMETHOD_REVIEW.md.

- На164строках максимальная разницаUDD−CF: **0,048585процентного пункта**.
- При предположении nearest6decimal rounding каждого исходногоmx±0.5e−6
  максимальная ширина propagated envelope — **0,000972процентного пункта**.
- Максимальное отличие произведения сохранённых авторскихqx от recomputedCF:
  **0,000235процентного пункта**. Это отдельный diagnostic, не officialUNqx.
- Проверены72age/sex chaining rows×2methods=**144равенства**
  Q(a,10)=Q(a,5)+(1−Q(a,5))Q(a+5,5),a40…75;164horizon comparisons10≥5.
  Возрастное снижение primarycurve не обнаружено; smoothing/исправление ставок
  не выполнялись и монотонность возраста не была условием допуска.

Rounding envelope **не CI**. Он не включает исходную неопределённостьWPP,
регистрации, unit bridge, calendar/treatment changes, heterogeneity или
индивидуальную calibration. РазницаCF/UDD тоже не является таким интервалом.
Малый эффект округления не доказывает высокую точность прогноза. Внутренние
85digit числа сохранены для арифметического replay, не как85значная точность
эпидемиологических данных.

## Источник и границы проверки

Source — уже сохранённый `aging_biomarkers@31fff9249237a33d9a10ac107d6710a04bfec517`,
`data/mortality/mx_qx_usa_deu_rus.csv`, SHA256
`59104bce947c3ab1132893e71decf946af3725c151fa5ff6f64bef440d4b54a2`.
Все4pinned source/context files проверяются **до CSVparse**. Выбрано100строк:
RUS/2023/F,M/годовые возраста40…89, seriesestimate. Нет person-level записей.

Зафиксированная provenance — UNWPP2024/mx1dt через redistribution package;
локальные документы задают annual deaths/person-year. В39 originalRDA checksum,
его pinned upstream version и конкретное numeric scaling package→CSV **не
аутентифицированы заново**. Поэтому это явно названный local-derived reference
при заявленных единицах, не новая сверенная officialUN life table. Official UN
methodology подтверждает определение central rate, а не автоматически bytes
этого CSV. Сохранённыйqx рассчитан авторами из originalmx отдельно и округлён
до6decimals; это не импортированная officialUNq.

[UNWPP2024 методика](https://population.un.org/wpp/assets/Files/WPP2024_Methodology-Report_Final.pdf)
и [country sources](https://population.un.org/wpp/assets/Files/WPP2024_Data_Sources.pdf)
подтверждают estimate/model происхождение. ДляRUS2023 нельзя утверждать
«зарегистрированные смерти2023»: перечисленные российские mortality inputs
заканчиваются раньше и включают корректировки. Latest2023 означает только
последний год в **пиннутом локальном CSV**, не globallylatestRussianrelease.
Методический источникUDD: [ONS national life-table guide](https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/methodologies/guidetocalculatingnationallifetables).
Подробности/точные страницы/ограничения —METHOD_REVIEW.md.

Общая популяция уже содержит людей с заболеваниями, лечение и все причины
смерти. Этот reference нельзя прибавлять к postdiagnosis all-cause mortality
из38 или использовать как disease-free hazard. Для общей calibration biomarker
модели нужны совместимые representative riskfactor/state distribution и
независимые целевые исходы. В39 они не получены; AUC/Brier/calibration/utility
не оценивались. Панель биомаркеров не подбиралась по этим возрастным таблицам.

## Параллельная квалификация SHFM

Официальный updated2013 Seattle Heart Failure Model рассмотрен отдельно:
[SHFM_REVIEW.md](SHFM_REVIEW.md), [rights note](SHFM_RIGHTS.md).
Подтверждена опубликованная однолетняя связь survival/score с0.9604. Численные
baseline2–5лет, полный обновлённый predictor/therapy/units/reference contract не
получены. Однолетняя точка не определяет пятилетнюю кривую; нельзя смешивать
старую exponential extrapolation и updated2013baseline. SHFM не реализована.

Read-only executorGET
`https://depts.washington.edu/shfm/update.php` для staticHTML/source checksum
получил **denied/stopped** `<urlopen error Tunnel connection failed:403 Forbidden>`
на proxyCONNECT. Предшествовавший webtextreview успешен; GET не повторялся через
headers/mirrors/transport. Role selector, acceptance, calculator requests и
контакты не выполнялись. Open software licence не установлена; публичная
читаемость не даёт прав на копирование/перевод исходного calculatorcode.
Missing source contract — scientific prerequisite pending, **не pendingapproval**.

Полный upstreamWPP rawRDA в39 не запрашивался; metadata navigation restrictions
описаны вMETHOD_REVIEW.md и не обходились. Отсутствие rawchecksum не названо
отказом в загрузке originaldataset, которая не выполнялась. Все источники/
артефакты предыдущих итераций остаются неизменными; reportv3/презентация не
редактировались. Новая ветка локальная, без push/PR/merge/deploy.

## Воспроизведение и статусы проверок

Из корня checkout, primaryCLI толькоstdlib:

```bash
python -m unittest tests.test_period39 tests.test_analyze39 -q
python -m research.mortality.period39.analyze --out /tmp/period39_NEW_run
MPLCONFIGDIR=/tmp/period39_mpl_cache XDG_CACHE_HOME=/tmp/period39_font_cache python -m research.mortality.period39.plot --results /tmp/period39_NEW_run/RESULTS_SUMMARY.json --svg /tmp/period39_NEW_reference.svg --png-preview /tmp/period39_NEW_reference.png
```

Нужны новые output paths и4локальных файла поSOURCE_MANIFEST.json; автоматических
загрузок нет. Chart dependency —Matplotlib3.10.8 в использованном окружении,
не обязательна для основного расчёта. Первые cachewarnings readonlyhome были
устранены заданием writable cache environment при повторном **только chart**
экспорте; SVG/PNG оказались побайтно идентичны. Численный164grid не повторялся.
ГрафикPNG визуально осмотрен, provenance вCHART_PROVENANCE.json.

Протокол: `c0a570e0ac40de94da8a05b53f5e9d80aaf7aacf`.
Scientific code: `151fe91b1540b91926baf3b3cc782cc7515ae55b`.
**27/27 unit tests PASS отдельноPython3.12.14 и3.13.5**; primary calculation и
chart export завершеныexit0. Code/metadata pins проверены до/после запуска;
4external source hashes —доparse и затем независимым current/Git аудитом.
RESULT_PROVENANCE.json содержит полные команды/hashes. GitHub CI для
этой локальной неопубликованной ветки не запускался.

[Независимый аудит](INDEPENDENT_AUDIT.md): **10295/10295PASS**. Проверены14code
pins,3неизменённых protocol files,4source files current/Git31fff,100source keys,
164rows/CSV и85digit вычисления,72×2chaining identities,164horizon comparisons,
20выбранных строк,metadata/scope и chart/sourceSHA linkage. Пять повреждённых
численных/scope fixtures отклонены даже при согласованно обновлённых result/CSV/
log hashes. Полный positive audit проверял также графики; они исключались только
из mutation fixtures, специально проверяющих арифметику и смысл полей.

Root дополнительно запустил независимый verifier на уже сохранённых flat-results
подPython3.13.5:10295/10295PASS, результат побайтно идентичен. Verifier не
импортирует и не перезапускает producer, не скачивает данные. Это проверка
арифметики/артефактов, не independent clinical validation. Он читает сохранённые
27-test logs; сами эти текстовые логи не аутентифицируют binary или remoteCI.

```bash
python research/mortality/period39/verify_results.py --repo "$PWD" --results research/mortality/period39 --out /tmp/period39_NEW_audit.json
```

`--results` допускает original aggregate/ и flat layout; разные nested/flat
выходы при одновременном наличии блокируются. Для воспроизведения используется
зафиксированная локальная source layout изSOURCE_MANIFEST.json. Источник нельзя
молча заменить другим release/путём с неподтверждёнными bytes.

Следующий содержательный шаг — квалификация полного источника для индивидуальных
условных эффектов/исходного состояния и представительной российской смеси
факторов риска. Этот периодный ориентир даёт воспроизводимый контроль масштаба,
но сам по себе не завершает перенос модели или клиническую пригодность панели.
