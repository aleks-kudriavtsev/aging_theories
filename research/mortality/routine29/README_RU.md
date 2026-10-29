# Итерация 29: белки сверх renal/routine core в повторно используемой NHANES III

Добавочная ценность reported Cys/B2M/BTP **сохранилась по Brier** сверх заранее
заданного compact routine core. Изменение AUC остаётся неопределённым.
Это новая exploratory-проверка на повторно используемой когорте, не независимая
внешняя или клиническая валидация. План7800458 записан до новых сравнений.

## Сопоставимый результат

Общий complete-case домен для ВСЕХ шести arms: **3 730 участников, 563 смерти
за восемь лет**,40–75 лет на интервью. WTCYPEX6,49страт/98PSU,KishESS1030,04;
раннего цензурирования к8годам нет. Взвешенная наблюдаемая смертность7,93616%.

| Заново обученная модель на этих же3730участниках | Brier8лет | AUC8лет |
|---|---:|---:|
| Simple core как28, новая оценка на общем домене29 |0.06567469|0.78521184|
| Simple core+3белка, тот же общий домен |0.06132253|0.81048400|
| Routine core |0.06369898|0.80073750|
| Routine core+3белка |0.06112195|0.81190804|

Основной заранее заданный контраст routine+3−routine:
ΔBrier **−0.0025770233**,95%paired full-nested-refit percentile interval
**−0.0036370857…−0.0008804442**;
ΔAUC **+0.0111705332**,interval **−0.0048276726…+0.0310748526**.
Все500парных повторов определены, ошибок основного nested refit0.
Выигрыш по правильному вероятностному score Brier поддерживается этой
внутренней оценкой; уверенное улучшение ранжирования по AUC не установлено.

Средний прогноз routine7,97660%,routine+3 7,99265%; bias+0,04045/+0,05649п.п.
OOF calibration intercept/slope:−0.12323/0.94030 и−0.16359/0.91894;
O/E0.99493/0.99293. Однозначного улучшения калибровки нет. Диагностические
fit не меняют оцениваемые прогнозы. Для каждой primary arm intercept/slope
определены498/500bootstrap-повторов: их bounds условны по определённым
повторам, не полные nominal95%intervals. Основные Brier/AUC определены500/500.

## Новый домен, пропуски и недопустимость cross-domain сравнения

В28 было3857/612; новый домен3730/563. Удалены127участников, среди них49смертей.
Это не случайная новая выборка; complete-case отбор может быть информативным.
Сравнивать сохранённые абсолютные Brier/AUC28 с29 для доказательства эффекта
нельзя. Simple модели в таблице заново обучены/оценены на том же домене29.

В old28-common причины исключений29:99только UACR,16только HbA1c,8толькоBMI,
4совместно BMI+UACR; сумма127. В более широком eligible protein-weight домене
3899 measurement-level пропуски: originalclinical7,proteins35,BMI14,HbA1c16,
eGFR0,UACR107. Они перекрываются; точные patterns в RESULTS_SUMMARY.json.
Взвешенное покрытие общего домена —97,2875% eligible protein-weight домена;
дополнительная коррекция весов за complete-case отбор не производилась.

## Routine core и первичные определения

Одиннадцать predictor terms: прежние age/age²,sex,current/former smoking,SBP;
плюс BMI,HbA1c,eGFRlow/high иlog2UACR. Это компактная спецификация, не вся
возможная routine/clinical панель. Ни CHD/angina proxy, ни новый diabetes-history
индикатор не добавлены. HbA1c используется непрерывно, без диагностического
порога. Креатинин иeGFR одновременно не включались.

CDC CEPnote подтверждает SCr_cal(mg/dL)=0.960*rawCEP−0.184,190сывороточных
образцов, сравнение с IDMS-traceable методом. Формула применяется ТОЛЬКО к
сывороточному CEP, никогда к urinaryURP. Из положительного SCr_cal по
NIDDK2021 CKD-EPI creatinine-only формуле рассчитан eGFRcr; race/Cys terms нет.
Это оценка, не measuredGFR; калибровка исторического Scr не доказывает точность
GFR у каждого пациента. Два фиксированных терма min(eGFR,60)/15 и
max(eGFR−60,0)/30 задают renal форму без поиска knots. Предположения и источники:
EGFR_RATIONALE.md и RENAL_ASSAY_REVIEW.md.

UACR=100*UBP(ug/mL)/URP(mg/dL). До моделей исключён URP<10:15eligibleзаписей,
все с опубликованной statistically-suspect replacement7.9. UBP0.4 — опубликованная
замена нижеLOD0.5, сохраняется как reported, без повторной импутации; она есть
у69общих участников. Не объявляется aboveLOQ. Другие отсутствующие UACR также
остаются missing; сывороточная калибровка не переносится на urine creatinine.

BMPBMI — released derived/edited BMI, иногда с официальными заменами массы/
роста. В общем домене weightflags0/3/4 имеют3613/80/37записей; heightflags
0/1/2/3 —3671/36/5/18. Все BMI не называются непосредственно измеренными.
HbA1c GHP% и GHPMETH1..4 сохранены; обе технологии описаны как DCCTaligned,
это не новое подтверждение современной NGSP-сертификации. Missingcodes
8888/88888/888888 обработаны ДО формул. SAS FORMAT .1/Z6.1 не является
implied-decimal INPUT: никакого ошибочного деления сырых значений на10 нет.

## Неопределённость и заранее заданные проверки

Основной анализ:5outer×3inner whole-PSU CV,train-only transforms/scaling,
ridgeλ1/10/100 по innerweightedBrier. Те же случаи/folds/weights на каждой
стороне контраста.500paired rescaledPSU повторов из eligible3899 с фиксированным
common mask заново выполняют ВСЮ preprocessing/tuning/outerfit процедуру.
Интервал условен по fixed partition,формам и выбранным гипотезам; не включает
assay/eGFRcalibration uncertainty, историю выбора моделей или внешнюю переносимость.

Пять дополнительных outcome-blind PSU-разбиений, без выбора лучшего:
ΔBrier **−0.00243442…−0.00213726**,отрицателен во всех5;
ΔAUC **+0.00798347…+0.01408097**,положителен во всех5.
Это5зависимых point-sensitivities, не5валидаций и не дополнительный CI.
Они не смешаны с bootstrap-интервалом основного fixed-разбиения.

Другие заранее помеченные exploratory варианты, тот же основной домен:
ранги трёх белков Brier0.06262762/AUC0.80028359;
training-only1/99%clipping Brier0.06117314/AUC0.80909307.
Ранговый вариант сохраняет point Brier-выигрыш, но не AUC-выигрыш относительно
routine. Варианты не выбраны по результату, отдельного superiority-вывода нет.

Metadata sensitivity: BMIflags0/0,GHPMETH1..3,UBP>0.4,поддомен3332/496смертей.
Обе модели заново оценены на нём: Brier0.06313351→0.06064449,
Δ−0.00248901; AUC0.79325886→0.80541627,Δ+0.01215741.
Это within-subdomain point-сравнение безbootstrap; абсолютные scores этого
поддомена не сравниваются с основным как доказательство улучшения качества.
Фильтр не является индивидуальной QC/LOQ-сертификацией.

## Проверяемость и CI

Получены14официальных источников. Новые LAB29314строк/6fields иEXAM31311/4
двумя reader согласовали175884и125244ячеек точно; прежниеclinical200500 также
согласованы. ASCIIдесятичные literals сохранены. Источники сохраняют разную
неиспользуемую trailing ширину по официальномуMISSOVER, ни одной записи не
пропущено. Потоковый sourceSHA совпадает с14-file manifest. EXAM195090426bytes
в пределах заранее обоснованного220MBcap; SAS никогда не исполняется.

53новых теста прошли локально иCI Python3.11/3.13. Научный codecommit
b8d05a0b2f6a534a98fa397a9f6862829d256df0:
https://github.com/aleks-kudriavtsev/aging_theories/actions/runs/40897873302
Empirical123531809031 иworkbench40897873199 success. Итоговые docs commits
пропускают CI, не повторяя модельные расчёты; это не новый CI на docs head.
PR121 browser40894133261 иmodel40894133044 проверены:success.

Независимый аудит **181/181**:80научных PythonSHA,план/metadata/source hashes,
90наборов fold-параметров(30main+50partitions+10metadata),арифметика/domain,
500CSVстрок. Основные интервалы независимо пересчитаны stdlib sorting+linear
percentiles, совпадение точное. Не выполнялись новые model fits,raw download
или ZIP. Дополнительный локальный GET только lab.sas metadata был denied
сетевым proxy403; повторной независимой проверки sourcebytes вне runner нет.
Источник-хеши сопоставлены между независимым stdout,manifest иreader audit;
это аудит агрегатов/провенанса, не отдельная переработка всех raw записей.

Сохранены RESULTS_SUMMARY.json,SOURCE_MANIFEST.json,FOLD_PARAMETERS.json,
CALIBRATION_BINS.json,primary_bootstrap.csv и INDEPENDENT_VERIFICATION.json/md.
В них нет participant ID,исходных записей или индивидуальных прогнозов.
ZIPartifact12908156978 доступен до12ноября2026; он не скачивался.
Все500paired deltas сохранены независимо от срока артефакта. Точныеrun/commit/
artifact данные —RESULT_PROVENANCE.json. Прежние binary28 и frozen models
проверены и сохранены; main,приложение,граф и coefficients не изменялись.

## Запуск

```bash
python -m pip install numpy==2.3.5 pandas==2.2.3 scipy==1.17.0 scikit-learn==1.8.0
python -m unittest discover -s tests -p 'test_routine29*.py' -v
python -m research.mortality.routine29.retrieve --out /tmp/new_routine29_inputs
OPENBLAS_NUM_THREADS=1 python -m research.mortality.routine29.analyze --source /tmp/new_routine29_inputs --out /tmp/new_routine29_results
python tools/verify_routine29_aggregates.py --job-log successful_job.log --out /tmp/new_routine29_audit
```

Дляauditнужен сохранённый plain-text stdout успешного job,неZIP. Все output
каталоги новые. Старые результаты не перезаписываются. Аккаунты,соглашения,
restricted data,контакты с держателями,merge/approvals/deployment не использовались.

Следующий обоснованный шаг — заранее заморозить полный новый preprocessing/
алгоритм и проверить на иной разрешённой когорте с теми же измерениями,
либо отдельно изучить недостающие routine/inflammation признаки. Повторный
подбор лучшего белка/разбиения на этих же OOF-прогнозах не станет независимой
проверкой. QC/LOQ,хранение,публичные perturbed followup времена и complete-case
selection по-прежнему ограничивают интерпретацию. Клиническая готовность
и доказанная биологическая полезность каждого белка не заявляются.
