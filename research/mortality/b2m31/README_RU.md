# Итерация 31: B2M сверх routine core на повторно используемом SSCARD_A

Добавление reported B2M к заранее выбранному routine core поддерживается
внутренней оценкой Brier и AUC. **Добавочный выигрыш B2M сверх routine + Cys
не установлен убедительно:** оба вторичных интервала пересекают ноль.
Это не выбор оптимальной панели, не confirmatory significance и не независимая
внешняя валидация frozen30. Выбор B2M направлен предыдущими результатами.

## Общий домен и сопоставимые результаты

Публичный исправленный SSCARD_A, биологический baseline NHANES 1999–2004,
не assay 2018–2020. Screening age 40–75. Из 6546 допустимых по linkage/весу
участников общий complete-case / measurable-flag00 домен — **6044 участника,
586 смертей за восемь лет**. Все четыре модели заново обучены и оценены на
одних людях, одинаковых фолдах и весах. Доли по циклам: 1889 /2131 /2024.
43 masked strata, 87 masked PSU; Kish effective n 3906.20. Минимальный
follow-up живых —177 месяцев, больше 96. Нулевые месяцы у умерших остаются
событиями. Взвешенная наблюдаемая смертность —7.449934%.

| Модель на одном общем домене | Brier, 8 лет | AUC, 8 лет |
|---|---:|---:|
| Routine core | 0.06069525 | 0.79106972 |
| Routine + B2M | 0.05928575 | 0.80262638 |
| Routine + Cys | 0.05965619 | 0.80218076 |
| Routine + Cys + B2M | 0.05922456 | 0.80445402 |

| Заданный контраст | ΔBrier | 95% refit percentile interval | ΔAUC | 95% refit percentile interval |
|---|---:|---|---:|---|
| Первичный: B2M сверх routine | -0.0014094962 | [-0.0022145495; -0.0006569203] | 0.0115566625 | [0.0043594763; 0.0235438811] |
| Вторичный: B2M сверх routine + Cys | -0.0004316326 | [-0.0010227332; 0.0001512836] | 0.0022732590 | [-0.0027857527; 0.0081626673] |

Все 500 парных повторов определены для каждого Brier/AUC контраста; ошибок
прогнозных refits во всех четырёх arms —0. Это conditional masked-PSU
refit intervals с фиксированными фолдами, формами и выбранной гипотезой,
не поправка за историю выбора маркера и не официальный NCHS variance estimator.

Нельзя заключать, что B2M лучше Cys, что два белка образуют оптимальную панель
или что отсутствие убедительного вторичного выигрыша делает B2M бесполезным
вообще. Сравнение B2M-vs-Cys и выбор победителя не являлись контрастами плана.

## Калибровка и неопределённые диагностические fits

Средний прогноз: routine7.461053%, routine+B2M7.467536%, routine+Cys7.475548%,
routine+Cys+B2M7.468930%; они не используются для рекалибровки. Уверенное
улучшение калибровки не установлено. Point diagnostic slopes и число
определённых bootstrap-диагностик:

| Модель | Slope | Intercept/slope valid | Offset valid |
|---|---:|---:|---:|
| Routine core | 0.980562 | 481/500 | 498/500 |
| Routine + B2M | 0.966047 | 489/500 | 498/500 |
| Routine + Cys | 0.972526 | 492/500 | 500/500 |
| Routine + Cys + B2M | 0.971128 | 487/500 | 497/500 |

Для intercept/slope это 19/11/8/13 undefined повторов; offset —2/2/0/3.
Такие интервалы условны по определённым fits, не полные nominal95% intervals.
Отсутствие ошибок основных прогнозных refits не означает, что все 500
калибровочных диагностик определены. Полные значения/статусы — в JSON и wideCSV.

## Веса, определения и missingness

Использованы пересозданные в июле 2022 SSCARD веса, взятые из самого SSCARD_A:
**(2/3)*WTSSCB4Y для 1999–2002; (1/3)*WTSSCB2Y для 2003–2004**. Не применялось
WTSSCB2Y/3 ко всем циклам, обычные MEC-веса или отдельная нормализация каждого
цикла перед scoring. Код/test/review подтверждают это правило; дополнительные
поцикловые raw/scaled weight sums не экспортированы и не были заново
вычислены независимым аудитором из XPT. Источники: SOURCE_WEIGHTS_REVIEW.md.

Уникальный masked cluster — (SDMVSTRA, SDMVPSU). Полный DEMO inventory
совпал с eligible inventory: 42 страты по 2PSU, одна (stratum1) с 3PSU.
Публичные inventories экспортированы как design metadata, не participant IDs.

Core имеет 11 фиксированных термов: age/age², sex, current/former smoking,
research mean первых трёх положительных SBP, BMI, HbA1c, два creatinine-only
CKD-EPI2021 eGFR терма и log2UACR. Креатинин 1999:1.013*rawLBXSCR+0.147;
2001LBDSCR/2003LBXSCR —без дополнительной коррекции. Ни формула CEP III,
ни старые cystatin calibration, ни Cys-eGFR здесь не используются. UACR:
100*URXUMA(ug/mL)/URXUCR(mg/dL), без III-specific urine floor10. Нет
диабетического/CHD proxy или одновременного raw creatinine и eGFR predictor.
Точные первичные определения — MEASUREMENT_REVIEW.md, CREATININE_REVIEW.md.

Из eligible домена исключены 502 человека; weighted coverage общего домена
92.940710%. Перекрывающиеся причины: missing routine474, missing/nonpositive
protein9, Cys not-flag0 14, B2M not-flag0 33. Их нельзя суммировать как
непересекающиеся причины; patterns сохранены в JSON. Missing routine включает
SBP286, BMI148, UACR84, creatinine/eGFR22, smoking11, HbA1c4. Коррекция весов
за этот complete-case/LOD отбор не выполнялась; исключённая популяция может
отличаться по риску. Death count среди исключённых отдельно не экспортирован.

## Assay/LOD sensitivity

Оба белка — reported serum mg/L, Siemens Dimension Vista1500. Primary требует
оба флага 0; это within-detection metadata, **не индивидуальный QC/LOQ pass**.
IFCC/ERM traceability SSCYST по release metadata не доказана. SSPRIS в общем
домене:6001 pristine/43non-pristine; они сохранены. Специальное предупреждение
для non-pristine glycated albumin не переносилось автоматически на эти белки.

Единственная заранее заданная sensitivity включает flags0/1 обоих белков,
сохраняя опубликованные LLOD/sqrt2 fills. Flag2 исключён из обоих доменов,
поскольку точная upper-LOD fill policy не описана. Sensitivity domain:
**6054 участника /586 смертей**, плюс 10 участников,0 событий. На его общем
домене primary ΔBrier−0.0014225866/ΔAUC+0.0117679503; secondary
ΔBrier−0.0004287565/ΔAUC+0.0021495655. Это point-only within-domain pairs,
без bootstrap CI; абсолютные scores разных доменов не сравниваются как
доказательство улучшения. Никаких дополнительных subgroup/transform поисков.

## Протокол, intake amendment и воспроизводимость

Первоначальный план commit **6774fbe600acad62cacbb9bece716f22ccef3560**
записан до сравнений. Первый scientific code db36608 остановился при intake:
50 файлов загрузились, два mortality PDF получили HTTP403; model step SKIPPED.
До любых сравнений зафиксирована узкая поправка **f27c7b47157db7fbaf31dc998bbb5c55c1c4f7e8**.
Она сохраняет уже проверенные официальные определения/CDC R specification и
три точных historical mortality SHA, не меняет outcome/weights/features/models
и не ослабляет readers. PDF не повторялись и не заменялись обходным URL.
Подробно: SOURCE_INTAKE_AMENDMENT.md, INTAKE_REVIEW.md, SOURCE_INTAKE_FAILURE.json.

Финальный run **40908831254**:66 tests на Python3.11/3.13 и empirical success;
workbench40908831274 success. Codehead f27c7b4 проверен через PR merge checkout
**a752d10a251f4bc0a381fe44d9d9c0a117d64f6e**; все 90 scientific hashes совпадают.
Последующий documentation-only head не считается новым scientific CI.
PR123 browser40904909281 и model40904909011 также завершились success.

Верифицированы **50 обязательных inputs**, все 25 исторических data SHA,
включая SSCARD_A actual official anchor URL; два ранее отклонённых PDF отдельно
отмечены denied_not_retried, без fake bytes/SHA. Текущее точное согласие двух
mortality readers:31126 строк/249008 ячеек. Положительность/схема/полные published
SSCARD marginals и corrected-weight counts проверены runtime guards. Сами
полные observed source marginals отдельно не экспортированы: независимый
аудит не заявляет второе чтение XPT или повторный расчёт индивидуальных весов.

Независимый численный аудит — **221/221**, включая все 40 интервальных полей
из 500 wideCSV строк, два narrow paired CSV,40 fold records и calibration bins,
полные bytes SOURCE_MANIFEST,90 scientific hashes и 84 сохранённых прежних.
Точные семь marked stdout blocks сохранены без JSON/CSV перезаписи.
Индивидуальные записи/ID/прогнозы не экспортируются; ZIP не скачивался.
Source/weight/domain аудит — **398 проверок,0 отказов**; его границы описаны
в SOURCE_INDEPENDENT_VERIFICATION.json/md.

```bash
python -m pip install numpy==2.3.5 pandas==2.2.3 scipy==1.17.0 scikit-learn==1.8.0
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_b2m31*.py' -v
python -m research.mortality.b2m31.retrieve --out /tmp/new_b2m31_inputs
OPENBLAS_NUM_THREADS=1 python -m research.mortality.b2m31.analyze --source /tmp/new_b2m31_inputs --out /tmp/new_b2m31_results
python tools/verify_b2m31_aggregates.py --log /tmp/saved_b2m31_job.log --ci-jobs /tmp/b2m31_independent_ci_jobs.json --ci-evidence-dir /tmp --out /tmp/new_b2m31_audit
python tools/audit_b2m31_sources.py --repo . --log /tmp/saved_b2m31_job.log --reference research/mortality/b2m31/SOURCE_EXPECTED_MANIFEST_REFERENCE.json --out /tmp/new_b2m31_source_audit
```

Для аудита нужны также три completed-job plaintext logs, названные
b2m31_independent_job_JOBID.log, в указанном ci-evidence-dir. Это разрешённые
GitHub log/metadata reads, не raw downloads. Все analysis output каталоги новые;
старые результаты не перезаписываются. Полный manifest/CSV доступны в repo,
поэтому результаты не зависят только от срока временного artifact. Точные
SHA/run/job/artifact expiration — RESULT_PROVENANCE.json.

## Границы результата и следующий шаг

SSCARD_A ранее импортирован в fe414dc1957d8c6041dce48e77cf58ee65b30d28;
Cys fitted в M4_cystatinC, B2M ранее loaded. Файлы/участники не новые для
программы; пересечение с прежним **fitted** M4 domain отдельно не восстановлено,
и нельзя называть оставшихся untouched holdout. Старые age40–79/10year/кардиальные
assay домены, веса и модельные формы отличаются; их числа здесь не сравниваются.
Это новый исследовательский B2M-контраст в повторной когорте, а не перенос
коэффициентов frozen30 или новая full3 когорта: BTP отсутствует.

Неопределённость не включает assay standardization/storage, предшествующий
выбор маркера, complete-case selection или внешнюю переносимость; pooled
whole-PSU CV не является forward calendar validation. Клиническая готовность,
аналитическая валидация assays и causal biology не установлены. Неопределённый
вторичный результат явно отражён; отсутствие убедительного выигрыша сверх Cys
не доказывает общую бесполезность B2M. Из него не выбирается новый победитель панели.

Следующий приоритет — независимый разрешённый cohort/assay bridge по freeze30
qualification, который остаётся отдельным unresolved blocker. Без нового
доступа допустима лишь новая заранее описанная sensitivity/reproducibility
проверка; дальнейший подбор маркеров на этих же OOF outcomes не станет внешней
валидацией. Frozen30/предыдущие модели/main/UI сохранены. Аккаунтов, соглашений,
restricted requests, контактов, merge/approvals/deployment не было.
