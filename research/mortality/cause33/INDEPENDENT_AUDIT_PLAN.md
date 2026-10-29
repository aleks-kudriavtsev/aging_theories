# Cause33A: независимый план аудита count/support gate

Подготовлен2026-10-29 до новых cause counts, по immutable gate protocol `82eacda12eced3c0105019fee09236b647821f45`. Репозиторий не меняется. Никаких raw retrieval, XPT/ASCII/ZIP/PDF GET, model fits, Fine–Gray или stageB comparisons. Future gate JSON/schema будет разобран после завершения официального source-qualified run; старые или незавершённые jobs не опрашиваются.

## Зафиксированные ожидания

- Пул1999/2001/2003,8лет/96месяцев inclusive отMEC, общий inherited flag00 routine11+Cys/B2M domain. Известные6044участника/586allcause deaths — только frozen domain consistency checks.
- ПуловыеW6:1999/2001→(2/3)WTSSCB4Y,2003→(1/3)WTSSCB2Y; native temporal weights32 здесь не являются score weights.
- Все10public underlying UCOD_LEADING groups001..010 иunknown отдельно; contributing flags не определяютcause. Для этих baseline лет доступны10groups; поздние2015–2018baseline cycles не поддерживаются.
- Allcause at/before96months: disjoint known10+unknown сумма raw events и weighted event totals. Unknown early cause блокирует любые later comparisons; unknown после96месяцев не делает8y target неизвестным, но reported separately. Alive-with-cause/нецелый илинеизвестныйcode/неверныйtime/status/schema — structural error.
- Дляcause k иfold cell: **N=event_k+known_noncase_k+unknown_early**. Known noncases включают competing known early deaths, deathsafter96 иobserved8y survivors; unknown early никогдане кодируется negative. ПриU=0 этообычноеN−event_k.
- Rank001..009 by descending weighted8y event total, затем descendingraw8y event count, затемascendingcode. Пять первыхопределяютсяДОsupport;010 остаётсявaccounting,неранжируется. Top5−unsupported неbackfill lower ranks. Возможны0..5eligible.
- Outer5salt`b2m31:outer`; inner3within eachoutertrain salt`b2m31:inner:{outer}`. Hash mapping восстанавливаетсяstdSHA256 на full eligible87maskedPSUs доcommonfilter. Full=eligible=common inventory,43strata,one3PSUstratum сохраняется.

## Арифметика, которую аудитор пересчитает независимо

1. Exactstdout blocks→savedJSONbytes+SHA, plan/LOCK/endpointmetadata, all95oldscientific hashes+полныйactualnewset; no mutation frozen30/old32. CI commit/merge-ref scope иsynthetictestcount — только final completed run.
2. Sourcewrapper↔50retrieved input fingerprints↔25oldpublicdataSHA↔retainedreadercontract;2oldPDFdenials separatelydenied_not_retried,безhash/placeholder/retry. Это reconciliation exported runtime hashes,неповторноелокальноеrawrehash.
3. Known10+unknown overall/percycle accounting: события,weightedtotals,общиеdenominators; percycle totals→pooled. All-followup cause counts≥ соответствующих8y count; earlyunknown≤all-timeunknown; deathsafter96 иalive учтеныотдельно.
4. Weighted incidence=sum_event_w/common_sum_w; share=sum_event_w/allcause_event_sum_w. Всеcauses имеютодинаковыйdomain denominator. Event Kish=(sum_event_w)^2/sum_event_w2; concentration=max_event_w/sum_event_w. Нулевыеdenominators→честноеundefined/None. Еслиmaxweightнеэкспортирован,fractionможнопроверитьлишьнаграницы/код,неclaimполноенезависимоевосстановление.
5. Независимаясортировка9cause records,exact ranks1..9,top5fixed,eligibility=priority∧support∧globalgate; no010/no backfill.
6. Длякаждойcause40fold cells:5outertrain+5outereval+15innertrain+15innereval. Outer eval totals суммируютсявfull; eachoutertrain+eval=full. Within eachoutertrain inner eval3суммируютсявoutertrain; eachinnertrain+eval=outertrain. Проверятьevent/noncase/unknown/n отдельно, анепредполагатьunknown0.
7. 83elementary supportchecks/cause:3global(eventrows≥50,knownnoncaserows≥50,event-bearingPSUs≥20),20outerthresholdchecks(5eval×event/noncase≥5,5train×event/noncase≥20),60innerchecks(15eval×event/noncase≥3,15train×event/noncase≥10). Каждаяfailure должнаиметьточныйcontext/observed/minimum; no post-count weakenedgate.
8. StageA flags: models_fitted=False,predictive_metrics_computed=False,no λ-tuning/model bootstrap/prediction arrays,no summingindependentbinaryprobabilities,FG/joint deferred. Counts/weightedburden/support неявляютсяmarkerbenefit,causality,power/EPV/calibrationcertificate илиexternalvalidation.

## Пределы независимости

Аудитор использует stdlib и агрегаты, без import scientific gate/pipeline длявычислений. Hash maps/ranking/support/ratio math восстанавливаютсяотдельно. Нетповторногоподсчётаперсональныхcauses изraw records; исправностьисточникаосновываетсянаretainedprimarydefinitions,sourcepins иruntime readers. Редкиеcells/flags сохраняютсякакaggregates, individual IDs/labels/predictions неэкспортируются.

Подготовительные материалы: `/tmp/cause33A_audit_expectations.md` ототдельногоагента; будущийstandaloneauditor `/tmp/verify33A.py`. Actualsource/codecommit/jobIDs иschema adapter будутзакрепленыпослеготовностиhelper иfinalCI,безновыхfitsилиданных.
