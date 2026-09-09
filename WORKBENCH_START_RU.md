# Mortality Workbench 0.7 — локальный запуск

Исследовательский интерфейс к зафиксированной модели общей смертности. Не предназначен для диагностики, выбора лечения или оценки современного индивидуального риска. Инфекции и внешние причины включены. Поддержка: исторические данные США, возраст 40–79 лет, горизонты 1/5/10 лет. Расчёт для RU/DE и 20 лет заблокирован.

## Без установки

Из корня распакованного архива или репозитория:

```bash
python3 -m research.mortality.workbench doctor
python3 -m research.mortality.workbench serve
```

Открыть `http://127.0.0.1:8765`, нажать «Заполнить синтетический пример», затем «Проверить и рассчитать форму». Сервер слушает только loopback. Остановка — Ctrl+C. В Windows вместо `python3` может потребоваться `py -3`.

## Установка Python-пакета

В отдельной виртуальной среде:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
mortality-workbench doctor
mortality-workbench serve
```

В Windows активировать `.venv\Scripts\activate`. Сборка из исходников устанавливает инструменты сборки из `pyproject.toml`; самому приложению сторонние библиотеки не нужны. Для полностью автономной установки подготовленного wheel:

```bash
python -m pip install --no-index --no-deps mortality_workbench-0.7.0-py3-none-any.whl
```

## Файлы и пакетный расчёт

```bash
mortality-workbench demo --out reports/new-demo
mortality-workbench example --model M4_cystatinC > example.json
mortality-workbench predict example.json --out reports/new-expanded
mortality-workbench batch inputs.jsonl --out reports/new-batch
mortality-workbench evidence --out reports/new-evidence
```

Все выходные каталоги должны быть новыми. JSON-контракт остаётся версии `0.6.0`, версия приложения — `0.7.0`. Не подменять BNP/NT-proBNP, тропонин I/T и RDW-SD/CV. Пропущенные обязательные данные, неизвестные поля и несовместимые единицы блокируют расчёт. Проверка материала/метода означает проверку заявленного ввода, не независимую верификацию лаборатории.

GET `/api/health` возвращает `ready_research_only` или HTTP503. Это техническая готовность, не клиническая. После изменения ввода старый результат становится недействительным; ответ на ранее отправленный запрос не должен подменять новый контекст.

## Проверки и воспроизведение

```bash
python -S -m unittest discover -s tests -p 'test_workbench*.py' -v
python -m research.mortality.verify_workbench06 --sources NHANES_SOURCE --out parity.json
```

Первая команда — 88 автономных тестов приложения. Вторая требует pandas/numpy/scipy/sklearn и прежних официальных файлов NHANES; она сравнивает реализацию, а не создаёт новую независимую валидацию. Полный научный набор требует WHO, NHANES и результатов итерации04; без них пропуски нельзя выдавать за полный прогон.

Исходные данные участников не входят в приложение. Сервер не сохраняет запросы и не отправляет измерения во внешние сервисы. Созданные пользователем отчёты содержат его ввод: не публиковать их в GitHub. Встроенный `http.server` — локальный исследовательский сервер, не готовый публичный медицинский сервис.
