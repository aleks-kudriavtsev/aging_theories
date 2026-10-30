# Iteration39: period mortality reference — source and method review

Дата проверки:2026-10-30. Задача — описательный age/sex reference для России, исходный точный возраст40–80, основной горизонт5лет, дополнительный10лет. Это **до вычисления** новых многолетних вероятностей. В данном аудите прочитаны локальные extractor/docs, существующий source gate36 и публичные первичные методические материалы. CSV заново не разбирался: проверен только SHA256. Нет новых risk products, численных возрастных прогнозов, fits, GET наборов данных/ZIP, аккаунтов, соглашений, contacts или изменений репозитория.

## Решение и предел квалификации

Сохранённый CSV можно воспроизводимо использовать как **локальный производный агрегатный источник с заявленной annual central-rate единицей** и явно названным преобразованием после фиксации протокола39. Это не эквивалентно независимому восстановлению официальных UN life-table probabilities: rawRDA bytes, его фиксированная версия и численный unit bridge именно к исходному `mx1dt` не аутентифицированы. Новая проверка официальной методики подтверждает значение central mortality rate, но не заменяет отсутствующую побайтовую проверку upstream release.

Предлагаемая первичная трактовка результата: **period2023 synthetic-cohort reference under a stated within-age mortality assumption**. Нельзя называть его наблюдённым пятилетним исходом российской когорты, прогнозом2023→2028/2033, персональной клинической вероятностью, baseline hazard для disease-free человека или новой российской калибровкой biomarker-модели. Статус series=`estimate` означает оценку за исторический период, не индивидуальное наблюдение.

## 1. Проверенные локальные данные и provenance

Репозиторий `/workspace/aging_biomarkers` имеет HEAD`31fff9249237a33d9a10ac107d6710a04bfec517`; working tree при проверке чист. Файл `data/mortality/mx_qx_usa_deu_rus.csv` имеет ожидаемый SHA`59104bce947c3ab1132893e71decf946af3725c151fa5ff6f64bef440d4b54a2`.

Схема, уже проверенная source gate36: `country,sex,year,age,mx,qx,series`. Зафиксированный домен RUS/2023/F,M/возраст40…89 содержит100 строк, без missing/duplicate/unexpected keys; finite положительные mx, qx∈[0,1), series=`estimate`. Эти утверждения прочитаны из сохранённого `SOURCE_GATE_RESULTS.json`, не пересчитаны по CSV в39. Пол B и person-year/population/death числители отсутствуют. Простое среднее F/M не является общеполовой вероятностью.

Фактический extractor:

- источник — mutable `https://raw.githubusercontent.com/PPgp/wpp2024/master/data/mx1dt.rda`, без upstream commit/SHA;
- `mx=float(r[mxM|mxF])`; дополнительного деления на1000 нет;
- сохраняет `mx=round(original_mx,6)`;
- **отдельно из original_mx**, а не из уже округлённого mx, считает `qx=round(1-exp(-original_mx),6)`;
- по умолчанию оставляет годы≤2023 и только M/F. Не был запущен в этом аудите.

Следствия: storedqx — авторская производная колонка, не импортированная официальная UNqx; исходное неокруглённое mx из неё не восстанавливается однозначно. Точное равенство `storedqx == 1-exp(-storedmx)` не требуется. Прежний gate допускает1e−6 для двух округлений; это tolerance проверки provenance, не confidence interval. Для нового референса метод должен явно потреблять mx и самостоятельно преобразовывать его по фиксированному правилу, сохраняя старый CSV/результаты неизменными. Старый qx допустим как provenance cross-check, а не автоматически другой независимый источник.

## 2. Что подтверждает первичная UN-документация

[UN WPP2024 Methodology Report](https://population.un.org/wpp/assets/Files/WPP2024_Methodology-Report_Final.pdf), July2024, UN DESA/POP/2024/DC/NO.10: разделI.C описывает разные способы оценки mortality/life tables по качеству исходных данных. AnnexA, printedp65/PDFp73, описывает single-age `1mx` через deaths с offset person-years of exposure; это central death rate, не само по себе event probability. Годовая сетка age/time1×1 введена в WPP2022. UN-методика не устанавливает для нашего производного CSV, что `1-exp(-mx)` является официальной WPPqx, и не подтверждает numerical scale преобразования package→localCSV.

[UN WPP2024 Data Sources](https://population.un.org/wpp/assets/Files/WPP2024_Data_Sources.pdf), August2024, UN DESA/POP/2024/DC/NO.11, **Russian Federation printedp205/PDFp214**: mortality inputs включают official estimates through2017, registered deaths by age/sex through2020, adjustments older-age mortality и international estimates through2021; adult mortality также корректируется на under-registration. Population estimates through2023 находятся в **другом** абзаце. Поэтому наличие country/year2023 не означает непосредственно зарегистрированную смертность2023. Это UN assessment/estimate, не новый national-register extract.

[UNData WPP2024 source metadata](https://data.un.org/Data.aspx?d=PopDiv&f=variableID%3A77%3BcrID%3A116&q=cambodia) разделяет estimates1950–2023 и projections2024–2100. Страница использована только как официальное описание revision/timeline; никакие значения чужой страны не использованы.

[Первичный README авторов package](https://github.com/PPgp/wpp2024) подтверждает `mx1dt` long format и колонки mxM/mxF/mxB с historical estimates/projections в одном файле. Он не даёт достаточного численного unit bridge к локальным bytes. Это документация redistribution package, а не новый официальный RDA release audit.

## 3. m, q и сила смертности — разные величины

Официальные определения и формула для взрослых подтверждаются [ONS Guide to calculating national life tables](https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/methodologies/guidetocalculatingnationallifetables),2019, и [NCHS U.S. Decennial Life Tables for1999–2001: Methodology of the United States Life Tables](https://www.cdc.gov/nchs/data/nvsr/nvsr57/nvsr57_04.pdf),2008-08-27. ONS приводит `q=2m/(2+m)` для возрастов>1. NCHS явно связывает `m=d/L` с uniform distribution of deaths для этой формулы. Это методические примеры life-table construction; они не являются американскими/британскими данными для России.

[Statistics Canada methods2023](https://www150.statcan.gc.ca/n1/pub/84-538-x/84-538-x2023001-eng.htm), Step3/6, независимо подтверждает actuarial formula и separation factor0.5 для взрослых. Их dimensionless factor не следует смешивать с нашим `a_x`, выраженным в годах.

Ниже — самостоятельный алгебраический вывод, а не расчёт по российским исходам.

Для интервала `[x,x+n)`:

```text
nd_x = l_x - l_(x+n)       — life-table deaths
nq_x = nd_x/l_x             — вероятность смерти в интервале
nL_x                       — person-years, прожитые в интервале
nm_x = nd_x/nL_x            — central death rate; единица1/year
na_x                       — среднее число лет, прожитых в интервале
                              теми, кто умер в нём
```

Если `na_x` выражено в годах, то из
`nL_x = n*l_(x+n) + na_x*nd_x` получается

```text
nq_x = n*nm_x / [1 + (n-na_x)*nm_x].
```

Для n=1 и uniform distribution of deaths (UDD), `a_x=0.5`:

```text
q_UDD(m) = m/(1+0.5*m) = 2*m/(2+m).
```

Отдельное допущение — постоянная сила `mu` внутри года возраста. Тогда
`l_(x+t)=l_x*exp(-mu*t)`, deaths=`l_x*(1-exp(-mu))`, exposure=`l_x*(1-exp(-mu))/mu`, поэтому центральный `m=mu` и

```text
q_CF(m) = 1-exp(-m).
```

Постоянная сила **не означает** `a_x=0.5`: её implied conditional mean time-to-death равен `1/m - 1/(exp(m)-1)` при m>0 и лишь стремится к0.5 при m→0. Комментарий `mx_to_qx()` об “implicit ax=0.5” математически неверен. Старый файл не менялся; новое исследование должно явно исправить трактовку в собственном протоколе/отчёте.

Оба правила — разные within-interval модели. Neither m nor q is automatically a patient-specific hazard. Нельзя утверждать, что UDD — установленная истинная поправка к данным UN для каждого возраста40–89. При отсутствии официальных age-specific ax в CSV его применение является предзаданной sensitivity-моделью. Для UDD необходимо проверять допустимость q∈[0,1]; большое m нельзя молча clipping-исправлять. CF вычислять численно устойчиво как `-expm1(-m)`, а survival/log-survival сохранять согласованно.

## 4. Period versus cohort и точные возрастные границы

Периодная таблица применяет age-specific rates одного периода к гипотетическому поколению. Когортная таблица следовала бы одному поколению через разные календарные годы. Это официальное различие изложено во введении NCHS2008, printedp1, и Statistics Canada2023. Фиксация year2023 на всех достигнутых возрастах создаёт period reference; она не предполагает, что реально наблюдавшаяся российская когорта состарилась с40 до90 в2023.

Индексная конвенция для протокола39:

- старт — **точный день рождения** целого возраста a; не произвольное число дней после последнего дня рождения;
- горизонт h соответствует полуоткрытому возрастному интервалу `[a,a+h)`;
- потребляются h полных интервалов возраста: a,a+1,…,a+h−1, всегда year2023;
- a=40,h=5 требует40…44; возраст45 — конец окна, а не дополнительный год экспозиции;
- a=80,h=10 требует80…89 и заканчивается в90; mx90 для этого окна не требуется;
- исходный диапазон40–80 включительно содержит старт80. Для fractional age или unknown birthday необходим отдельный протокол; молчаливое округление/интерполяция не разрешены этим дизайном.

Формула только для последующего выполнения после freeze:

```text
period_survival(a,h,2023) = product_{j=0..h-1}(1-q_(a+j,2023)).
```

В этом аудите она не вычислялась. Символический индекс подчёркивает отличие от будущего calendar-cohort маршрута с меняющимся year; последние rates отсутствуют в frozen historicalCSV.

## 5. Шесть знаков: sensitivity, а не статистический CI

`round(original_mx,6)` даёт представительное значение mhat и интервал квантования примерно

```text
delta = 0.5e-6
m_low = max(0, mhat-delta)
m_high = mhat+delta.
```

Закрытые границы — консервативная конвенция для округления к ближайшему; binary-floating/tie behaviour допускает только отдельно названный машинный tolerance, не расширение ради прохождения данных. Это последняя экспортная точность, не полная upstream precision/error model. Отдельно округлённый qx не нужно выдавать за истинный q или сужать им интервал без протокола совместного rounding-consistency.

Оба допустимых преобразования монотонны по m; после freeze компонентные нижние/верхние m можно протянуть через выбранное правило и многолетнюю формулу. Одинаковый знак крайних perturbations по всем возрастам даёт консервативный **rounding envelope**, не случайные независимые draw и не95%CI. В39 пока не получены ни эти численные envelopes, ни age/sex predictions.

Не смешивать три разные вещи:

1. CF versus UDD — чувствительность к within-age shape assumption;
2. mhat±0.5e−6 — чувствительность к экспортному округлению;
3. UN-estimation/revision uncertainty, death under-registration, age misstatement, treatment/calendar changes и patient heterogeneity — этой sensitivity не охвачены. В CSV нет полных intervals/covariance для этого слоя.

## 6. Рекомендуемые fail-closed проверки до будущего выполнения

1. Сверить локальный SOURCE_MANIFEST и все source/code/plan hashes **до parsing**; не запускать downloader или старые builders.
2. Повторить точную схему и RUS2023/F,M/40…89 inventory gate; не подменять missing нулём, чужим годом или другим полом. Фильтр baseline age40…80 и required attained-age support40…89 — разные операции.
3. Числа конечны, не bool; год/возраст имеют точную целочисленную семантику; series=`estimate`. Единица mx зафиксирована как locally documented deaths/person-year; upstream-authenticated flag остаётся false.
4. Основной CF и вторичный UDD именованы заранее; storedqx не выбирает между ними. Явно блокировать недопустимое q/потерянный возраст/ошибочную единицу, не clip и не substitute.
5. Вывести только агрегированный descriptive reference и его named sensitivity, со scope flags: `individual_prediction_validated=false`, `disease_free_baseline_identified=false`, `cohort_forecast=false`, `original_RDA_hash_verified=false`. Не комбинировать этот слой с заболеванием/стадиями или markerHR в39 без нового предзаданного identification/transport шага.

Это рекомендации к ещё не исполненному протоколу, не утверждение об уже выполненной реализации39.

## 7. Manifest-ready local context hashes

Все SHA256 ниже получены чтением bytes, без парсинга данных и без модификаций.

| Absolute path | SHA256 |
|---|---|
| `/workspace/aging_biomarkers/data/mortality/mx_qx_usa_deu_rus.csv` | `59104bce947c3ab1132893e71decf946af3725c151fa5ff6f64bef440d4b54a2` |
| `/workspace/aging_biomarkers/tools/fetch_wpp_mortality.py` | `2b652175319688d7f0932db44d2a732ba5f349e897c08727604cab06379a8519` |
| `/workspace/aging_biomarkers/data/mortality/README.md` | `f54f34b69e67e3c25ce0a8b4cbf95b6640d42a2d3dd0d4abbbf2bb72b90101ed` |
| `/workspace/aging_biomarkers/docs/12_cross_country_mortality.md` | `ba87e883a6a149c49d1c74d00d0eb2f105c2e3c46e290c603368fed761a55268` |
| `/workspace/aging_theories_registry37_local/research/mortality/modular36/SOURCE_MANIFEST.json` | `94cee585675447cce935fd88452a48192a7c047683cd8c2bdb56c3c5e78d2625` |
| `/workspace/aging_theories_registry37_local/research/mortality/modular36/SOURCE_GATE_RESULTS.json` | `f67ec3f08fa0e93e03c7fff1b1291e0f6df206451b4501dddaa0f0832e2f8608` |
| `/workspace/aging_theories_registry37_local/research/mortality/modular36/qualify_sources.py` | `491bcc1dbc1eac3a2912bf1a221151e4aa85d7cf8628ef3aaf82b137f5ded217` |
| `/workspace/aging_theories_registry37_local/research/mortality/modular36/RUSSIAN_SOURCE_AUDIT.md` | `028b9eb532fe3e95828291165f08580cdd2ab5aa9508e9ddcad72ca24d8bcafa` |

`SOURCE_GATE_RESULTS.json` is prior36 output, not a newly generated39 result. Репликацию100 qualified rows здесь не объявлять новой независимой проверкой.

## 8. Доступ: точные ограничения

- Прочитаны официальные UN methodological/data-source PDF и primary metadata, не новые datasets. Полные upstream data bytes не запрашивались.
- Предположенный старый путь UN glossary и `master/man/mx.Rd` дали404; это не source evidence. ПутьCRAN HTML сообщил inaccessible. GitHub package rootREADME читался; клики в `tree/main/man` и `tree/main/data-raw` вернули **restricted URL**, после чего этот путь остановлен. API/raw/mirror обходов не было.
- Поэтому завершено то, что доступно: единицы/методологические определения и локальная provenance/rounding логика. Независимая численная unit/byte equivalence rawRDA→CSV остаётся **not verified**, не “denied original dataset”: dataset retrieval не предпринимался и не был разрешён этой задачей.

Новых российских вероятностей, combined risks, ранжирования, clinical validation или обучения не выполнено.
