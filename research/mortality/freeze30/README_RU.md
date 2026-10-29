# Итерация 30: замороженный исследовательский алгоритм и готовность к переносу

Созданы два новых final-training артефакта: **routine core** и **routine core +
три белка**. Коэффициенты, порядок признаков, единицы, преобразования,
нормализация, подбор ridge и endpoint полностью определены. Научный CI успешен.
**Новая независимая внешняя валидация не выполнена:** ни один из 12 проверенных
кандидатов/релизов не прошёл все условия доступа и научной совместимости.

Это отдельная восьмилетняя бинарная ridge-logistic модель. Она не заменяет
прежний cause-specific engine, frozen models, main, UI или калькулятор.
Снижение Brier в исследовании 29 остаётся прежним внутренним результатом;
финальное обучение 30 не создаёт новой оценки AUC, Brier или калибровки.

## Что именно заморожено

План 5cbf90e записан до финального обучения. Его конкретизация источников
SOURCE_POLICY.md и исполняемая реализация 16e2eb8 также предшествуют обучению.
LOCK.json связывает план, CONTRACT.json, источник и результаты 29; перед чтением
данных проверяются все 80 исторических научных Python-хешей.

Training domain совпадает с 29: **3730 участников /563 смерти к 96 месяцам**,
возраст 40–75 на интервью, первый exam NHANES III, WTCYPEX6,49 страт/98PSU.
Complete-case отбор прежний; новое расширение популяции не заявляется.
Признаки: age/age², sex, current/former smoking, SBP, BMI, HbA1c,
два фиксированных eGFRcr terms и log2UACR. Во второй модели добавлены log2
original reported Cys/B2M/BTP. Никакого выбора отдельного белка, rank/clipping
варианта или лучшего разбиения в 30 нет.

Три whole-PSU inner folds, salt `freeze30:final-inner:v1`, lambda1/10/100;
взвешенный inner Brier используется только для настройки обучения, при разнице
не более 1e-12 выбирается большая lambda. Все преобразования обучаются внутри
inner training. Затем preprocessing и коэффициенты заново оцениваются на всём
training domain. Выбраны **lambda10 для routine и 1 для routine_all3**.
Финальные максимальные градиенты:7.40404e-7 и 2.63044e-8. Веса likelihood
нормируются к n обучения; intercept не штрафуется.

Два запуска на одинаковых входах и коде дали **побайтово одинаковые**
FROZEN_MODELS.json и FREEZE_MANIFEST.json. Внутри runner adapter воспроизвёл
численные прогнозы существующих primitives на 3730 случаях для каждой модели:
максимальное абсолютное расхождение **0**. Прогнозы остались временными
внутренними массивами; сохранены только параметры и агрегированные проверки.
Inner Brier и objective в артефакте — диагностические величины обучения,
не validation scores. Внешние performance metrics и risk CI не вычислялись.

## Строгий adapter и границы применимости

CONTRACT.json задаёт точные поля/единицы/категории и единственный mapping
`NHANESIII_FIRST_EXAM_REPORTED_V1`. Adapter принимает только полные записи
одного общего домена для обеих моделей. Неверные/лишние поля, неизвестный
mapping, единицы, null/nonfinite/bool numbers, sentinel codes, неподдержанные
возраст/метаданные или неположительные измерения отклоняются.

Исторический raw serum CEP преобразуется ровно один раз:
`SCr_cal=0.960*CEP−0.184 mg/dL`, затем creatinine-only CKD-EPI2021.
Современный IDMS-креатинин нельзя выдавать за raw CEP; сывороточная коррекция
никогда не применяется к моче. UACR=100*UBP(ug/mL)/URP(mg/dL), URP>=10.
Released UBP0.4 и официальные BMI edits сохранены с их исходными значениями.
Они не становятся точными above-LOQ измерениями или индивидуальным QC_PASS.
Современный стандартизованный Cys, NPX/RFU и повторный визит не преобразуются
автоматически в исторические reported concentrations.

`load_artifact` проверяет JSON, manifest и byte hashes; `score` возвращает
вероятности только в памяти. **Scoring CLI и автоматического экспорта записей
нет.** Этот API не разрешает публикацию индивидуальных прогнозов. Успешная
проверка схемы не подтверждает происхождение данных, права доступа или
взаимозаменяемость анализов. Hash не защищает от замены одновременно файла
и его manifest; проверяемый reference задаёт repository commit.

## Конкретная qualification, без имитации внешнего анализа

Проверка началась с актуальных external18/hrs24/cohort27 roadmap и веток.
Полные первичные ссылки, определения, пределы проверки и статусы действий:
QUALIFICATION_COHORTS_RU.md, QUALIFICATION_PUBLIC_RU.md, SSCARD_HISTORY_RU.md.
Машинные metadata contracts: CANDIDATE_MANIFEST.json; их исполняемые blockers:
QUALIFICATION_RESULTS.json. Неизвестные поля намеренно не заполнены догадками.

| Путь | Решающий blocker сейчас |
|---|---|
| ARIC visit4,1996–98 | Все 3 serum assays/UACR известны; same-visit HbA1c не подтверждён. BioLINCC approval/RMDA не получены. Современный assay требует bridge. |
| MESA-Kidney,2012–14 | Малый ancillary domain, не весь MESA. Пакет events2020 не обеспечивает 8 лет для поздних визитов; same-visit core и ancillary access не квалифицированы. |
| CHS / AGES-Kidney / CRIC | Контролируемый либо не подтверждённый доступ; exact joint panel, baseline и endpoint остаются частично неизвестны. |
| HRS2016VBS | Person-level AI-анализ запрещён политикой провайдера; полный 8-летний endpoint и вся панель не подтверждены. |
| UK Biobank | Новые заявки приостановлены; существующий допуск не предъявлен. Все 3 targets есть в NPX, это не serum mg/L. |
| SSCARD_A1999–2004 | Открытый релиз Cys+B2M, BTP отсутствует. Источник и участники уже использовались; Cys уже fitted M4, B2M ранее загружен. |
| NHANESIII first/repeat /2009–10 | Первый визит использован; второй — те же люди; для 2009–10 full3 в проверенных каталогах не найден. |
| UCI CKD | Открытый dataset не содержит нужные белки или mortality endpoint. |

Реестр **approved external mappings пуст**. Даже полностью заполненные
metadata-декларации не включают scoring. Для будущего переноса нужны точный
допущенный release, подтверждённый same-visit domain, assay bridge, отдельный
reviewed mapping commit и outcome protocol до любых внешних сравнений.
Нет новой внешней выборки, event count или validation metrics30.

## Проверки и воспроизведение

60 новых тестов прошли локально и на Python3.11/3.13. Научный codecommit
`16e2eb8fe0abd51b351a6c3775a7180ebdd25354`, run40903512074: synthetic и
final-training success; workbench40903512040 success. Научный CI относится
к этому codecommit, не к последующему documentation-only head. PR122
browser-regressions40900913543 и model-check40900913507 — success.

13 scientific/doc source files совпали с 29 по SHA точно. Live index.htm
изменился, но все resolved release URLs остались прежними; оба indexSHA
сохранены. Разрешение этого варианта каталога было зафиксировано до обучения
в SOURCE_POLICY.md. Изменение measurement/SAS/codebook байтов не допускалось.

```bash
python -m pip install numpy==2.3.5 pandas==2.2.3 scipy==1.17.0 scikit-learn==1.8.0
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_freeze30*.py' -v
python tools/check_freeze30_qualification.py
python -m research.mortality.routine29.retrieve --out /tmp/new_freeze30_inputs
OPENBLAS_NUM_THREADS=1 python -m research.mortality.freeze30.freeze --source /tmp/new_freeze30_inputs --out /tmp/new_freeze30_results
```

Каталоги должны быть новыми. Копии параметров/manifest сохранены из двух
отмеченных JSON-блоков успешного authorized job stdout; ZIP не запрашивался.
Artifact API metadata и точные fileSHA — в RESULT_PROVENANCE.json.
Независимый аудит: **106/106 проверок**,84 научных Python-хеша, сохранность
80 исторических хешей, оба вывода artifact/manifest и все12 qualification
результатов. Независимые формулы на шести вымышленных профилях совпали со
scorer до1.11e-16. Полный предел проверки — INDEPENDENT_VERIFICATION.json/md.
Новые байты retrieval manifest из runner не экспортированы: его SHA сверялся
между aggregate payloads, а не повторным локальным хешированием того файла.
Actual source bytes также не загружались независимым аудитором. Это не вторая
переработка raw cohort. Детерминизм подтверждён в одной pinned среде, без
гарантии побайтовой идентичности на любых платформах/версиях.

Для повторения аудита нужны четыре сохранённых authorized CI metadata/log
файла: `freeze30_independent_job.log`, `freeze30_jobs.json`,
`freeze30_synthetic311.log`, `freeze30_synthetic313.log`. ZIP не нужен:

```bash
python tools/verify_freeze30_aggregates.py --evidence-dir /tmp/saved_freeze30_logs --out /tmp/new_freeze30_audit
```

## Ограничения и следующий разрешённый шаг

Сохраняются reported-value QC/LOQ противоречия, неизвестная индивидуальная
пригодность образцов, хранение, historical assay/calibration uncertainty и
информативный complete-case отбор. eGFR — оценка, не measured GFR. Нет новых
clinical-readiness или external-transport claims. Прежняя неопределённость 29
не становится CI индивидуального риска frozen coefficients30.

Новая независимая эмпирика **заблокирована именно отсутствием квалифицированной
разрешённой когорты**, а не сбоем final training. Для ARIC нужны допустимые
человек/учреждение/страна, реальное одобрение/release и решение same-visit HbA1c;
для HRS AI-microdata путь остаётся запрещённым даже при login. Никаких заявок,
аккаунтов, соглашений или контактов не инициировано; это not_requested /
not_granted, не pending approval. Ошибки metadata403/reader-size отдельно
описаны в qualification: их не следует выдавать за denial доступа к person data.

Без нового доступа возможен только отдельный заранее заданный B2M-контраст
на уже известном SSCARD_A, с правильными 2022 весами, assay/LOD правилами и
новым shared domain. Это исследовательская работа на повторной когорте,
не проверка трёхбелкового frozen artifact и не независимый holdout. Для
внешней проверки полного алгоритма сначала требуется разрешить реальные
access/measurement/endpoint blockers. Слияний, approvals и deployment нет.
