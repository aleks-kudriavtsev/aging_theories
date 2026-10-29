# Forward32: независимый аудит источников, весов и домена

**Пройден: 380 проверок, отказов 0.** Проверены агрегаты, код и экспортированные манифесты; записи участников, raw-source bytes, ZIP и сетевые источники не открывались, refits не выполнялись.

Scientific head `59b0626767c0169107652cc7b7d51a416c8e7668`; checkout журнала — PR merge `def39d99ab86bb6f7b022554bdb9f7cb077e3d27`. Совпали 14 LOCK metadata SHA, 90 historical Python SHA и все 95 exported scientific SHA. Четыре JSON совпадают с marked job exports; шесть вспомогательных файлов совпадают с RESULT output SHA. Run/job указаны по переданному журналу и информации root; отдельного сетевого CI lookup не было.

Ровно 50 обязательных source entries; все 25 data SHA совпали с прежними pins по именам. Заявленные байты всех источников: 54,415,075; данных: 50,591,670. Свежий base manifest и wrapper32 связаны своими хешами, а старый manifest31 отдельно сохранён в LOCK. Два прежних PDF с HTTP403 отмечены denied_not_retried; новых запросов/хешей этих PDF не заявлено.

SSCARD runtime href: `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.xpt`. SHA `63c8342280dee3bd67c1322141ef42bc6a350e479cbaeb69fc4474a67b14796d`; размер 2,279,120 байтов. Полные observed marginals (оба flags0/1/2/missing, pristine и оба веса) точно совпали с заранее закреплёнными таблицами July2022 для 21890 строк. Это независимая сверка экспортированных source aggregates, а не повторное чтение raw records. Serum, mg/L и Dimension Vista1500 сохранены; LOD не назван LOQ.

Native train weight — WTSSCB4Y, late — WTSSCB2Y. Для каждого eligible/common периода и цикла проверены source/native суммы, originalW6 factors, sum of squared weights, Kish effective n и coverage. Rowwise max difference 0 — экспортированный diagnostic. Late W4 полностью missing; его сумма0 отражает sum(skipna), а не индивидуальные нулевые веса. Estimands условны на полученный complete-case measurable domain; национальная W6 estimate не заявлена.

| Часть / домен | n | События / несобытия | Сумма native weights | Kish n | Event-bearing PSU |
|---|---:|---:|---:|---:|---:|
| train:eligible | 4323 | 460 / 3863 | 105612581.312926 | 2755.537552 | 57 |
| train:common | 4020 | 383 / 3637 | 99016388.114702 | 2572.544991 | 57 |
| late:eligible | 2223 | 241 / 1982 | 111565329.830487 | 1464.898524 | 30 |
| late:common | 2024 | 203 / 1821 | 101970999.131736 | 1333.692684 | 30 |

Coverage common/eligible: train 0.937543490404, late 0.914002578459. Primary flags00 сохранены во всех четырёх arms: train — 4020/383, late — 2024/203. LOD01 sensitivity итерации31 здесь не анализировалась.

Full = eligible = common inventory в каждой части: train — 28 strata / 57 PSU, late — 15 strata / 30 PSU. Единственная m3 — train stratum1=[1,2,3]. Participant overlap и masked-PSU overlap — 0. Точная early fold mapping независимо восстановлена через SHA256: по 19 PSU в каждом из трёх folds. Агрегаты inner scoring складываются в полный train domain; late в них не входит.

Gate min20 пройден для positive-weight event и nonevent **rows** в каждой части. Это не порог 20 PSU и не доказательство мощности или калибровочной точности. Реальные common event-bearing PSU: 57 train / 30 late; nonevent-bearing PSU также 57/30. Минимальный survivor follow-up — 201/177 месяцев, выше 96.

Ограничения: нет batch/individual assay date/storage duration; assay 2018–2020 не равен baseline 1999–2004. Cohorts и поздняя часть уже использовались, маркер выбран после прежних результатов. Это retrospective baseline-calendar stress test, не untouched holdout, не deployment as-of2003 и не независимая валидация frozen30. Все пределы record-level независимости перечислены в companion JSON. Численные интервалы bootstrap проверяет отдельный reviewer.
