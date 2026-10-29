# Итерация31: SSCARD_A и routine core — проверка измерений до протокола

Проверено **2026-10-29** по первичным официальным документам CDC/FDA и существующей истории программы. Только метаданные: raw/XPT, индивидуальные записи и новые сравнительные исходы не открывались; модели не обучались; файлы репозитория не изменялись. Подробная независимая проверка сывороточного креатинина: `/tmp/creatinine31.md`.

**Вывод:** допустим новый, явно названный исследовательский анализ **reported CysC+B2M** на повторно используемой когорте NHANES1999–2004. Цикл-специфическая обработка serum creatinine разрешает обоснованный creatinine-only eGFR core. Нельзя объявлять это переносом frozen30, клинической assay-валидацией, полной трёхбелковой панелью или независимым новым holdout. Для SSCARD_A не установлена lot/release-specific IFCC-трассируемость cystatin; она не выводится из года измерения или названия Siemens.

## Белки: значения, флаги и разумный primary domain

[Официальный SSCARD_A codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.htm): биологический serum baseline1999–2004, анализы2018–2020, Siemens DimensionVista1500 immunoassays. Релиз пересмотренJuly2022: ошибочно включённые записи удалены, сохранённые значения не менялись, созданы новые веса.

| Поле SSCARD_A | Определение / шкала |
|---|---|
| `SSCYST` | CystatinC, serum mg/L; опубликованные пределы0.23–8.00mg/L |
| `SSB2M` | Beta2-microglobulin, serum mg/L; пределы0.72–23.0mg/L |
| `SSCYSTL`,`SSB2ML` | 0=в пределах detection;1=нижеLLOD;2=вышеULOD;dot=missing |
| Значение при flag1 | Уже подставлено LLOD/√2; повторно не импутировать |
| `SSPRIS` | 1=pristine,0=non-pristine; не индивидуальный QC-pass |
| `WTSSCB2Y`,`WTSSCB4Y` | Специальные новые веса;4Y только1999–2002 |

BTP отсутствует в полном перечне SSCARD_A. Предупреждение об ухудшении non-pristine performance в codebook относится к glycated albumin; его нельзя автоматически распространять на оба белка или считать pristine достаточным подтверждением их пригодности.

**Рекомендация до новых результатов:** primary common domain требует оба белка positive finite и **оба флага0**, полный routine core, положительный правильный вес, eligible/определённый8-летний outcome. Все сравниваемые arms пересчитываются на тех же участниках. Primary transforms — `log2(SSCYST)` и `log2(SSB2M)` как released concentrations; никакой автокалибровки к SSCYSTAT/SSCYPC, клиппинга по частотным диапазонам или выбора белка по результатам.

Если нужен один ограниченный LOD sensitivity, **заранее** определить альтернативный domain: флаги0/1, positive finite released values, flag2 и unknown по-прежнему исключены. Flag1 трактуется как опубликованная подстановка, а не точная концентрация; обе модели заново оцениваются внутри этого домена. Не сравнивать абсолютные Brier/AUC разных доменов как эффект маркера. Возможную pristine-only проверку надо отдельно зафиксировать до результатов либо оставить агрегатным описанием; не добавлять её из-за неудобного primary результата.

Flags0 означают within-detection, **не** автоматически измерение выше валидированного LOQ и не клиническую пригодность каждого образца. При необычном сочетании flag0 и значения вне опубликованных assay limits сначала разбирать источник/округление, а не незаметно удалять запись. Primary domain определяет условную популяцию и может исключать крайние почечные состояния; нужны overlapping reason counts и изменение weighted coverage.

## Почему нельзя обещать современную стандартизацию Cys

SSCARD_A не сообщает точный reagent/calibrator lot, номер пересмотренной методики или ERM-DA471/IFCC-трассируемость этого анализа. Исторический [FDA K063351, DimensionVistaCYSC](https://www.accessdata.fda.gov/cdrh_docs/reviews/K063351.pdf), pp0–3,6, описывает particle-enhanced immunonephelometry, те же0.23–8.0mg/L и CYSC-калибратор с highly purified proteins. Это подтверждает существование иной исторической калибровки, **но не доказывает**, какой вариант был использован в SSCARD_A2018–20. ERM-DA470 в общей calibrator/control документации нельзя принять за ERM-DA471/IFCC cystatin standard.

Следовательно защищённая формулировка — **reported SiemensVista concentrations**. Ни коэффициент1.11, ни старые NHANES cystatin recalibration equations не применять без отдельного release-specific основания. Не рассчитывать eGFRcys и не использовать absolute clinical cystatin thresholds в этой итерации. Для log2-маркера в новом исследовательском fitted contrast это ограничение интерпретации и transport, а не основание выдавать численную доступность за QC. Межплатформенная сопоставимость с III остаётся неподтверждённой.

## Точная карта routine/source по циклам

`SEQN` — ключ one-to-one join, а `SDDSRVYR` должен согласоваться с циклом:1/2/3. Текущие codebook columns необходимо сверить с фактическим release schema при штатном импорте; здесь raw columns не просматривались.

| Компонент |1999–2000|2001–2002|2003–2004| Поля, единицы, определение |
|---|---|---|---|---|
| Demographics/design |[DEMO](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/DEMO.htm)|[DEMO_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/DEMO_B.htm)|[DEMO_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DEMO_C.htm)|`SEQN`,`SDDSRVYR`,`RIDSTATR`,`RIDAGEYR`,`RIAGENDR`,`SDMVSTRA`,`SDMVPSU` |
| Smoking |[SMQ](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SMQ.htm)|[SMQ_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/SMQ_B.htm)|[SMQ_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/SMQ_C.htm)|`SMQ020`,`SMQ040`; одинаковые ключевые коды |
| Blood pressure |[BPX](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/BPX.htm)|[BPX_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/BPX_B.htm)|[BPX_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/BPX_C.htm)|`BPXSY1`,`BPXSY2`,`BPXSY3`,дополнительно`BPXSY4`,mmHg |
| Body measures |[BMX](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/BMX.htm)|[BMX_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/BMX_B.htm)|[BMX_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/BMX_C.htm)|`BMXBMI`,kg/m²; released BMI из weight/standing-height |
| HbA1c |[LAB10](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/LAB10.htm)|[L10_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L10_B.htm)|[L10_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L10_C.htm)|`LBXGH`,%; whole blood, Primus boronate-affinity HPLC |
| Serum creatinine |[LAB18](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/LAB18.htm)|[L40_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L40_B.htm)|[L40_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L40_C.htm)|**`LBXSCR` / `LBDSCR` / `LBXSCR`**,mg/dL; см. следующую таблицу |
| Urine albumin+creatinine |[LAB16](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/LAB16.htm)|[L16_B](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L16_B.htm)|[L16_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L16_C.htm)|`URXUMA`,µg/mL;`URXUCR`,mg/dL; исходные, не SI aliases |
| Protein shared release |SSCARD_A|SSCARD_A|SSCARD_A|Один совместный1999–2004файл; не копировать/повторно объединять его строки как три независимых релиза |

`RIDAGEYR` — **age at household screening, adjudicated**, не NHANESIII `HSAGEIR` и не возраст на MEC. Для нового core можно заранее использовать screening age40–75 и тот же age в eGFR, но нужно так и назвать. `RIDAGEEX` — месяцы на examination; его нельзя незаметно подставить вместо `RIDAGEYR`. `RIAGENDR`:1=male,2=female; неизвестное не кодировать какfemale.

Smoking: `SMQ020=2`→never; при `SMQ020=1`, `SMQ040∈{1,2}`→current (every/some days), `SMQ040=3`→former. Refused7,don'tknow9,missing не превращать в never. Пропущенный по маршруту SMQ040 у SMQ020=2 не делает never автоматически missing. Указать именно cigarette history≥100 lifetime, а не всякий tobacco exposure.

Для SBP есть3, иногда4 попытки. `BPXSAR` в1999–02 — результат, сообщённый участнику с отдельным алгоритмом исключения первого, а не универсальный traditional mean. Защищённый компактный вариант для31 — заранее заданный research mean positive finite `BPXSY1..3`, минимум1; это соответствует прежней программе, но **не** объявляется эквивалентом III released `HAZMNK1R`. Если используется4-я попытка или исключение1-й, это отдельный fixed rule до результатов. Нельзя переключить алгоритм по Brier. Primary BMI — released positive finite `BMXBMI`; не переносить III `BMPWTFLG/BMPHTFLG` codes как флаги этого release.

`LBXGH` — released standardized HbA1c%, не glycated albumin и не IFCC mmol/mol. CDC2012notes прямо требуют original HbA1c data без withdrawn cross-over regression. Эти циклы имеют Primus method; III `GHPMETH` как обязательного поля здесь нет. Не добавлять diabetes-history дубликат вместо непрерывного HbA1c.

## Serum creatinine: строго разные правила1999/2001/2003

| Цикл | Raw mg/dL source | Дополнительная CDC calibration к IDMS-traceable reference |
|---|---|---|
|1999–2000|`LAB18.LBXSCR`|**`Scr=1.013*LBXSCR+0.147`** |
|2001–2002|`L40_B.LBDSCR`|**identity**: дополнительная поправка не нужна |
|2003–2004|`L40_C.LBXSCR`|**identity**: дополнительная поправка не нужна |

Направление1999 equation — исходный NHANES serum result→CCF Roche enzymatic, IDMS/NISTSRM967 traceability. Результат1mg/dL становится1.160mg/dL. В2001–02 уже выполнена межлабораторная harmonization в release; не применять её повторно. Сравнительные диагностические регрессии2001/2003 из analytic notes не являются указанием дополнительно исправить значения. **NHANESIII `.960*x-.184` здесь неприменима.** Serum correction никогда не применяется к`URXUCR`.

В1999generic[methodPDF](https://wwwn.cdc.gov/nchs/data/nhanes/public/1999/labmethods/lab18_met_biochemistry_profile.pdf) поле названоLBDSCR, тогда как current cycle codebook указываетLBXSCR. Это документационная несогласованность: при импорте надо assert intended actual column и остановиться при расхождении, не делать silent alias/fallback. См. также[2001method](https://wwwn.cdc.gov/nchs/data/nhanes/public/2001/labmethods/l40_b_met_creatinine.pdf),[2003method](https://wwwn.cdc.gov/nchs/data/nhanes/public/2003/labmethods/l40_c_met_creatinine.pdf).

Creatinine-only CKD-EPI2021 допустимо назвать **estimated** GFR по этой CDC calibration guidance; race и Cys не входят. Дляfemale k=.7,alpha=−.241,factor1.012; дляmale k=.9,alpha=−.302; формула `142*min(Scr/k,1)^alpha*max(Scr/k,1)^(-1.2)*.9938^age*sex_factor`. Только renal terms `min(eGFR,60)/15`,`max(eGFR-60,0)/30`; не дублировать raw/logScr.

## UACR и missingness

Из исходных единиц всех3 urine releases: **`UACR(mg/g)=100*URXUMA(µg/mL)/URXUCR(mg/dL)`**. Не использовать serum creatinine в знаменателе и не применять serum calibration к urine. `URXUMASI`mg/L численно равноµg/mL, но`URXUCRSI`µmol/L — другая шкала; использование обеих SI колонок требует иной формулы, поэтому предпочтительны raw-поля выше.

Urine — spot specimen; creatinine Jaffé/CX3, albumin fluorescent immunoassay. Для primary reported UACR требуются оба positive finite; dot/NaN→missing, никаких нормальных значений вместо пропусков. **Не переносить** III URP7.9/порог10mg/dL: LAB16 допускает urineCr2mg/dL, L16_B5mg/dL в официальных частотных таблицах. Их диапазоны не являются клиническими cutoff. В этих codebook нет белковых SSCARD-подобных LOD comment fields; low-end urine values сохраняют reported-value ограничения, а не новый QC-pass.

В routine fields этих циклов используются numerical/dot missing codes из соответствующих codebook; III sentinels8888/88888/888888 не являются универсальными здесь. Reader обязан различать столбцы и значения фактического release; неизвестные типы/коды исследуются до моделирования. Полный comparison domain фиксируется до оценки и одинаков для всех arms.

## Веса, исходная точка и8лет

SSCARD требует именно его пересмотренные subsample weights. [Официальный NCHS weighting tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx), раздел combining1999–2002/1999–2004, требует использовать four-year weight первых двух циклов. Для объединённой six-year оценки: **`(2/3)*WTSSCB4Y` для1999–2002, `(1/3)*WTSSCB2Y` для2003–2004**. Для отдельногоtrain1999–02 иtest2003–04 — соответственноWTSSCB4Y иWTSSCB2Y; единая константа внутри подвыборки при normalized ridge/weighted mean не меняет оценку, но relative pooling weights важны. Не подставлять WTCYPEX6,ordinaryMECweights или naive WT2Y/3 всем. Whole-PSU/design handling требует отдельной фиксации protocol поSDMVSTRA/SDMVPSU.

[2019public LMF dictionary](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf): joinпоSEQN, ELIGSTAT1, MORTSTAT0/1, **PERMTH_EXM** отMEC доdeath/31Dec2019. Binary death≤96months допустим **после** проверки на exact common domain: определеныstatus/time; alive имеют≥96months; deadafter96 — survived8y,неevent. Если естьalive<96илинепонятное время — stop,неприсваивать0; IPCW потребует отдельного протокола. Календарноbaseline1999–2004до2019даёт>8лет, но это не заменяет проверку released times.

[Public-use description](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf) предупреждает о подстановке отдельных follow-up times/causes ради конфиденциальности; vital status не подменён. Assay dates2018–20 — даты анализа сохранённых baseline образцов, а не начало наблюдения. Согласно[актуальной NCHS странице](https://www.cdc.gov/nchs/linked-data/mortality-files/), открытая связка2019;2022restricted/RDC здесь не используется. Cause-specific/CHD/angina признаки не нужны и proxy не подставляются.

## Повторное использование и границы вывода

`freeze30/SSCARD_HISTORY_RU.md` фиксирует загрузку SSCARD_A и fitted Cys model в benchmark commit`fe414dc1957d8c6041dce48e77cf58ee65b30d28`; B2M тогда был загружен, хотя отдельный fitted contrast не найден в проверенной истории. Эти же календарные участники/исходы ранее исследованы. Поэтому31—новая заранее определённая исследовательская проверка добавочного B2M/двухбелковой панели на **reused cohort**, не independent untouched validation. Изменившаяся дата измерения stored sample, новыйконтраст или другойfoldsalt не создают новых людей. Full3/ BTP оценка этим не выполняется.

Оставшиеся ограничения: точные joint n/события/LOD-by-cycle/missingness/weights coverage ещё не вычислялись в этом review; individual QC, long-storage effects и lot-specific protein standardization не доказаны. Однако определения и released measurable flags позволяют ограниченный reported-value prognostic анализ при описанном pre-outcome плане. Любой source/schema конфликт должен быть разрешён до сравнительных вычислений, без ослабления reader validation.
