# External18: владелец данных выполняет замороженную модель локально

Первый исполнимый пакет внешней проверки RU/DE. Реальные иностранные когорты в этот выпуск не входят; полученных RU/DE метрик нет. Протокол требует совместного согласования, а не является уже одобренной заявкой.

## Структура

- PROTOCOL_RU.md — первичный вопрос/контраст и критерии выполнения.
- COHORTS.json, SOURCES.md — проверенные маршруты SHIP/KORA/ЭССЕ и текущие статусы доступа.
- runner.py — preflight до чтения участников; отдельный исследовательский запуск, метрики только общей смертности.
- precision.py — повторяемость/междневная/межзапусковая компоненты на парных контрольных уровнях.
- METHOD_COMPARISON_RU.md — предлагаемый лабораторный дизайн.
- demo.py — только явно синтетическая проверка программы.

Нужны numpy, pandas, scipy и scikit-learn из существующего научного окружения. Исходный Workbench и его country gates не меняются. Синтетические строки не считаются внешней валидацией.

```bash
python -m research.mortality.external18.demo --out /tmp/ext18_synthetic
python -m research.mortality.external18.runner preflight --config /tmp/ext18_synthetic/config.json --protocol /tmp/ext18_synthetic/protocol.txt --dictionary /tmp/ext18_synthetic/dictionary.json --out /tmp/ext18_synthetic/preflight.json
python -m research.mortality.external18.runner evaluate --config /tmp/ext18_synthetic/config.json --protocol /tmp/ext18_synthetic/protocol.txt --dictionary /tmp/ext18_synthetic/dictionary.json --data /tmp/ext18_synthetic/synthetic.jsonl --out /tmp/ext18_synthetic/result.json
python -m unittest discover -s tests -p test_external18.py -v
```

## Настоящее исследование

У держателя данных: утвердить протокол и доступ, сформировать словарь с методами/единицами, преобразовать свою выгрузку в канонический JSONL без имён и контактов. Реальные исходные коды переменных SHIP/KORA/ЭССЕ в этом пакете НЕ угаданы: их должен подтвердить держатель. Использовать data_origin=authorized_cohort, 1000 bootstrap-повторов, SHA документов/данных и положительные независимые подтверждения из ATTESTATIONS. Передать `--approval LOCAL_APPROVAL_FILE` обеим командам. Program preflight проверяет наличие/хеш/подтверждения, но не подлинность правового основания.

В JSONL ровно: record_id, age_years, sex, clinical, measurements, death, followup_years, weight, stratum, cluster. В clinical: smoking, bp_treatment, diabetes_history, cvd_history, cancer_history. В measurements: bmi, sbp, hba1c, creatinine, uacr, albumin. Неполные данные — null, не0. Death — true/false с установленным vital status; неизвестный исход надо сначала согласовать. Если причина смерти неизвестна, но факт известен, death=true. Единицы/метод/матрица заданы отдельно в словаре; код не принимает неизвестные преобразования. Для HbA1c IFCC→NGSP используется официальное аффинное уравнение; другие формулы не угадываются.

Только baseline-выживаемость без отложенного входа, 3года, возраст40–79 в поддержке когорты. Репрезентативность, информативное цензурирование и выбор полного профиля требуют содержательной оценки. При оценке причин следующей версией понадобятся исходные ICD10, точная семантика и отдельные тесты.

Выход содержит агрегаты без record_id, индивидуальных измерений/прогнозов, без разбивок малых групп. Он всё равно проходит disclosure review держателя: проверка минимума10 событий и известных выживших НЕ доказывает конфиденциальность или достаточную статистическую точность. Путь сохранения обязан находиться в разрешённой среде, не в публичном git. Не отправляйте данные или подписанные разрешения в GitHub Actions/Vercel. Автоматической сетевой передачи нет.
