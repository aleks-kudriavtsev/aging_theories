# Изолированная публикация пакета39

Разрешение пользователя «слей с веткой main» относится к текущему пакету39.
Ветка подготовлена от main@adc4d5dd80c27b2b3215147f6bad97445b72916c:
только period39, два его тестовых файла, offline CI. Итерации35–38 не добавлены.
Все14 scientific code/metadata pins, исходные результаты, графики и historical
provenance сохранены побайтно. Исторические commits в manifest остаются historical
идентификаторами; их ancestry не выдумывается в новой ветке.

Для автономной проверки сохранены четыре уже доступных публичных агрегатных/
контекстных файла source_snapshot/ с исходными hashes. Это не новые загрузки,
не ZIP и не обход прежних сетевых запретов. Person-level данных нет. Источник:
UN DESA Population Division, World Population Prospects2024; заявленная
source licence CC BY3.0IGO; rawRDA equivalence по-прежнему не проверена.

verify_publication.py сохраняет независимую арифметику исходного verifier,
но проверяет реальные published Git blobs и pinned original manifest вместо
недоступных historical ancestor objects. Он явно указывает этот предел.
verify_results.py и старые audit JSON не изменены. Новая проверка имеет собственный
счётчик; не подменяет historical10295-check report.

python -m unittest tests.test_period39 tests.test_analyze39 -v
python research/mortality/period39/verify_publication.py --repo . --results research/mortality/period39 --out /tmp/period39_publication_audit.json

Оригинальный producer требует прежнюю source layout /workspace/aging_biomarkers.
Offline CI не вызывает fetch/extractor/producer: заново считает независимые
Decimal/Fraction references из source_snapshot и проверяет все164cells/CSV.
Клиническая валидация/изменение UI/моделей не входят в публикацию.
