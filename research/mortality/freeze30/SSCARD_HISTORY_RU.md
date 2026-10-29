# SSCARD_A: дополнительная проверка происхождения

2026-10-29; только код, сохранённые агрегаты, manifest и уже открытые официальные HTML-метаданные. Никакой XPT-запрос, поиск новых raw-файлов или просмотр индивидуальных записей не выполнялся.

## Подтверждённое прежнее использование

Коммит `fe414dc1957d8c6041dce48e77cf58ee65b30d28`, 2026-09-10, `Add trained NHANES mortality benchmark with temporal validation, calibrated-risk audit, published coefficients and 95-test verification`:

- `research/mortality/nhanes_benchmark.py:103–109` импортирует SSCARD_A из локального `1999_SSCARD_A.xpt`, поля `SSCYST`, `SSCYSTL`, `SSB2M`; соединяет их по SEQN; Cys участвует в domain restriction.
- `PANELS` включает `M4_cystatinC`, где к прежнему core+NTproBNP+troponin добавляется `log_cystatin`; преобразование задаётся в том же модуле. `nhanes04/model_bundle.json` содержит обученную модель `M4_cystatinC`. Следовательно Cys уже анализировался здесь как предиктор.
- `SSB2M` читается/соединяется, но отсутствует в `PANELS`, `CONTINUOUS` и преобразованиях этого benchmark. Это подтверждает prior data availability/loading, но само по себе не fitted B2M analysis.
- `nhanes04/SOURCE_DATA_SHA256.txt:15` фиксирует `1999_SSCARD_A.xpt`: `63c8342280dee3bd67c1322141ef42bc6a350e479cbaeb69fc4474a67b14796d`. Это прежний manifest hash, не новая загрузка или новая проверка bytes.

Коммит `663eec48e7002c4f3ff9f0b0c065109ac6bf487d` добавляет/поддерживает тот же импорт в `nhanes_cohort06.py:29–35`. Обе версии используют уже известные 1999–2004 календарные когорты.

Проверенный `git log --all -G 'SSB2M|log_b2m|log_beta2|beta2_microglobulin' -- 'research/mortality/*.py'` на доступной истории показывает прежние benchmark/import изменения и последующие 27/28 III-изменения. Отдельного fitted B2M-контраста на continuous NHANES до III в этой проверке не найдено. Это ограниченный аудит доступной Python-history, не утверждение об отсутствии анализа в недоступных/удалённых архивах.

**Корректная формулировка возможной следующей работы:** заранее определённая исследовательская оценка добавочного B2M либо Cys+B2M-контраста на повторно используемых SSCARD_A/NHANES 1999–2004. Cys, исходный релиз и участники не новые; трибелковая модель не проверяется, BTP отсутствует. Нельзя представлять этот путь как обнаружение новой панели или новой независимой когорты.

## Точный direct XPT href

Уже открытый [официальный общий каталог Laboratory Data](https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory) содержит строку 1999–2004 SSCARD_A Data [XPT - 2.2 MB], отображённую веб-инструментом под ссылкой 562. [Самостоятельный codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.htm) указывает имя `SSCARD_A.xpt`.

Однако рендер веб-инструмента раскрывает текст и номер ссылки, а не её raw href. Переход по 562 запросил бы сам XPT и поэтому **не выполнялся**. Сохранённый raw HTML с href локально не найден. Python metadata route ранее получил proxy403 и не повторялся. Поэтому `direct_url=null`, `direct_link_resolved=false`; подтверждены release metadata и download-link presence, но не точная прямая ссылка/новая доступность bytes. Подстановка `.htm`→`.xpt` была бы выводом по шаблону, а не проверенным извлечением href.

Это `unresolved` metadata-reader blocker, а не отказ CDC в доступе к SSCARD_A и не неизвестный исход выполненной raw-загрузки: raw-запрос вообще не предпринимался.
