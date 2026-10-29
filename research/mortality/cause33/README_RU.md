# Итерация33A: какие отдельные причины поддерживают сравнение маркеров

**Данные поддерживают следующий сравнительный протокол только для двух из
пяти ведущих определённых причин:** злокачественных новообразований002 и болезней
сердца001. Это source/count qualification, а не результат оценки белков.
Ни одной cause-specific модели, Brier/AUC, predictive calibration или модельного
bootstrap в33A не рассчитывалось. Совместная интеграция и Fine–Gray отложены.
Результаты32 сохранены в draftPR24/evidencePR125; оба checks PR125 успешны.

## Один домен и один горизонт

Повторно используемый NHANES1999–2004/SSCARD_A, baseline MEC, horizon96months
включительно. Eligible6546; общий flag00 complete-case domain6044 /586 смертей.
Unknown causes к96месяцам0; известный восьмилетний статус у всех. Общий denominator
одинаков для всех причин. Сохраняются11routine features, оба serum reported
mg/L протеина Cys/B2M и прежние measurement/creatinine/QC/LOD ограничения.
Это не независимая валидация, не новый holdout и не перенос frozen30.

Pooled corrected six-year weights:2/3WTSSCB4Y для1999–2002,1/3WTSSCB2Y для2003–04.
Это не32's native temporal estimands. Common weight sum100001258.45371333,
weighted all-cause events7450027.753817, incidence7.449934%.
Взвешенная common/eligible coverage92.940710%; веса не устраняют complete-case,
linkage/assay/flag selection. Рейтинг относится к этому received domain,
не к России/миру или безусловной national population. Стратификация43/87PSU,
включая третью PSU stratum1; full/eligible/common inventories совпадают.

## Ранг по взвешенному mortality burden, а не качеству модели

|Код|Опубликованная underlying cause|Events8|Взвешенная8y incidence|Доля weighted8y смертей|Event PSU|Support|В следующий протокол|
|---|---|---:|---:|---:|---:|---|---|
|002|Злокачественные новообразования|176|2.383497%|31.9935%|75|pass|да|
|001|Болезни сердца|148|1.882176%|25.2643%|67|pass|да|
|003|Хронические болезни нижних дыхательных путей|33|0.429771%|5.7688%|25|fail (2)|нет|
|004|Несчастные случаи|25|0.422768%|5.6748%|21|fail (7)|нет|
|005|Цереброваскулярные болезни|25|0.255916%|3.4351%|20|fail (7)|нет|
|007|Сахарный диабет|23|0.203248%|2.7282%|19|fail (9)|нет|
|008|Грипп и пневмония|13|0.126711%|1.7008%|13|fail (30)|нет|
|009|Нефрит, нефротический синдром и нефроз|10|0.100496%|1.3490%|8|fail (33)|нет|
|006|Болезнь Альцгеймера|4|0.066604%|0.8940%|4|fail (42)|нет|
|010|Все остальные причины|129|1.578747%|21.1914%|60|pass|нет|

Последняя строка010 — неоднородный остаток, **вне клинического рейтинга**.
Её129смертей/21.1914% weighted mortality остаются в accounting и competing
noncases. Неизвестная причина не подменяется010. Все десять групп суммируются
в586 событий и7450027.753817 weighted events; это сумма observed disjoint counts,
не сумма независимо fitted вероятностей. Ранги определённых причин:
002,001,003,004,005,007,008,009,006. Top5 выбирался до support, без backfill.

Первые два проходят все83проверки на причину. CLRD003:33<50событий и
в outer-evaluation fold3 только4<5событий. Accidents004:25<50 и ещё6 недостаточных
training/evaluation cells (минимальный outer-evaluation2). Cerebrovascular005:
25<50 и ещё6 недостаточных cells (outer-evaluation minimum3,inner-training8).
Точные40cells и каждое нарушение сохранены. Остальные не заменяют эти группы.
Порог50events/20eventPSU и фиксированные fold minima — заранее выбранный
computational/reporting guard, не power/EPV/calibration certificate и не
утверждение о бесполезности белков при редких причинах. Для002 event-weightKish
121.894, maximum single-event weight1.621%; для00196.217 и2.108%.

## Что делает gate и чего ещё нет

Протокол82eacda12eced3c0105019fee09236b647821f45 зафиксирован до новых counts.
Девять именованных UCOD001..009 доступны в этих baseline years; ограничение
2015–18 до001/002/010 относится к baseline cycle, не году смерти. Primary
underlying cause не заменена contributing diabetes/hypertension flags.
Heart001 не вся cardiovascular mortality, CLRD003 не изолированныйJ44,
renal009 не изолированнаяCKD. Public2019 cause/time perturbation и assumed
linkage statuses остаются ограничениями; vital status не perturbated.

Неизвестная common ранняя причина блокирует все сравнения и не кодируется0;
смерть после96месяцев имеет известный noncase status на этом горизонте.
Известная конкурирующая смерть<=96 — noncase отдельной binary CIF, не censoring
и не утверждение о выживании. При alive followup<96 pipeline останавливается.
Отдельные fixed-horizon probabilities не cause-specific intensities и не
когерентный совместный набор; суммирование/новая joint система запрещены.

Следующий разрешённый шаг — отдельный committed сравнительный протокол для
002/001: routine±B2M±Cys, primary B2M сверх routine и secondary сверх routine+Cys,
одинаковый domain и5×3 whole-PSU training-only pipeline, paired full-refit
uncertainty, calibration/rare-undefined/multiplicity limits. В33A эффекты
маркеров ещё **не оценены**. Маркер не закрепляется исключительно за одной
причиной. Недостаток support для ещё трёх priority causes относится к текущему
8y дизайну; пороги/горизонт/split не менялись после counts. Fine–Gray/совместная
интеграция остаются отложенными, старый five-cause hazard engine сохранён.

## Проверки, происхождение и запуск

Научный code4ff18aa4115fa0b0d93a0addf9f5c1e0c86a8a83;
реальный PRmerge checkout53002102a1af5a5c23d7df127dd137dbb0f1856b.
CI40924184727: empirical123625177367 и28tests Python3.11/3.13 success;
workbench40924184882success. Все99scientific hashes равны codehead; прежние95
модулей и18locked metadata неизменны. Docs/results-only commit сskipCI не
выдаётся за новый scientific run. PRE_COUNT_FREEZE описывает небольшое
дополнительное static registry validation до первого empirical run.

50official inputs/25historical data SHA и два совпадающих mortality readers.
Два старых PDF403 не запрашивались снова; нет новых PDF checksum claims.
Current dictionary review использует official indexed text и retained mapping,
не заново полученные PDF bytes. ZIP не скачивался: два exact JSON stdout blocks
сохранены и побайтно совпали с независимым прочтением completed GitHub log.

Независимый агрегатный аудит373/373PASS: accounting, ratios/Kish/concentration,
percycle/cell complement identities, independent hash folds, все830support
comparisons/400cells, ranking/no-backfill и source/code provenance. Это не
повторное чтение raw participants и не независимая reanalysis их исходов.
Новых прогнозов, clinical validity или независимой cohort validation нет.

```bash
python -m pip install numpy==2.3.5 pandas==2.2.3 scipy==1.17.0 scikit-learn==1.8.0
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_cause33*.py' -v
python -m research.mortality.cause33.retrieve --out /tmp/new_cause33_sources
OPENBLAS_NUM_THREADS=1 python -m research.mortality.cause33.qualify --source /tmp/new_cause33_sources --out /tmp/new_cause33_gate
python tools/verify_cause33_aggregates.py --repo . --log /tmp/cause33_success_job.log --science-commit 4ff18aa4115fa0b0d93a0addf9f5c1e0c86a8a83 --run-id 40924184727 --job-id 123625177367 --ci-jobs /tmp/cause33A_independent_ci_jobs.json --ci-log-dir /tmp --tests-per-python 28 --out /tmp/new_cause33_audit
```

Для auditor нужны сохранённые completed-job metadata и plaintext logs с именами
cause33_independent_job_JOBID.log. Aggregate artifacts остаются в Git; expiry
GitHub archive не удаляет evidence. Только агрегаты, без participant IDs/labels/
predictions/raw files; accounts/agreements/restricted data/contacts/merge/
approval/deploy/scheduled tasks отсутствуют.
