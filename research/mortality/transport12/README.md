# Transport12 — воспроизведение научного анализа

Пять широких групп первоначальной причины смерти, исторические США, возраст 40–79 лет. Приложение Workbench по умолчанию не переключено. Строгий расчёт не заполняет отсутствующий CRP или другие обязательные измерения.

## Полученные результаты

Обучение: 15 280 полных профилей 1999–2010. Основная проверка: 2666 участников 2011–2012 и 130 смертей за пять лет. Дополнительная проверка: 3013 участников 2013–2014 и 113 смертей за четыре года. Это разные горизонты; смерти нельзя объединять в один пятилетний показатель.

Основной многокатегориальный Brier: 0,07104 → 0,06871, разность −0,00233 с условным интервалом −0,00409…−0,00069. Во втором цикле интервал разности включает ноль. Общая AUC повышается в обоих циклах; редкие респираторные и цереброваскулярные исходы имеют недостаточную точность проверки. Полный отчёт: `../ITERATION_12_RU.md`.

Анализ без RDW и лейкоцитов — отдельная проверка чувствительности, добавленная по документации анализаторов до просмотра новых индивидуальных исходов, а не замена первичной модели лучшим результатом.

## Неизменяемые контрольные точки GitHub

- План: `55279f23e10b9318ceddfa266d8d404e717006ee`, `ANALYSIS_LOCK.md`.
- Первичная модель зафиксирована: `70e8cf3505b8910d7b0689c56158b00d123c6823`, `MODEL_FREEZE.json`.
- Первичная модель SHA-256: `40069d1a98a29308bee26e5dac895388f4c9008b983bc834ad8ac0a57fbcd094`.
- Дополнение по метаданным: `d52e6b472a6c069b52a1cee23185641fae897afd`, `ASSAY_AMENDMENT.md`.
- Вторичная модель зафиксирована: `be9c8f18894ae93e6779c49fcc24779b6b130100`, `SENSITIVITY_FREEZE.json`.
- Вторичная модель SHA-256: `7c7dc918f042ff44f141a9ad3b49e6573db7f9732c5eb9d93e29b09b0686fa24`.

Полные тексты этих контрольных точек находятся в ветке `research/transport-no-crp-12` репозитория `aleks-kudriavtsev/aging_theories`. Протокол не выдаётся за внешнюю регистрацию исследования.

## Автономный пример, без данных участников

Из корня исходников:

```bash
python3 -S -m research.mortality.replay_transport12 --demo --out synthetic12.json
python3 -S -m unittest discover -s tests -p 'test_transport_runtime12.py' -v
```

Это синтетический пример расчёта, не медицинское заключение. Для настоящего исследовательского ввода нужны точные единицы и материал. Доступные горизонты — 1, 4 и 5 лет. Россия, Германия, 20 лет, точная ИБС/N18 и индивидуальные доверительные интервалы не поддерживаются.

## Точное воспроизведение оценки сохранённых моделей

Индивидуальные исходники не входят в научный архив. Их получают из официальных источников CDC с зафиксированными отпечатками. `expected_development_sources.json` содержит 87 отпечатков файлов разработки; `source_manifest.json` в сопровождающем архиве — 54 новых файла данных/документации. Исходные URL также записаны в workflow `transport12-*-inputs.yml` и предыдущих workflow проекта.

Установлены numpy, pandas, scipy, scikit-learn; для проверок независимого обучения требуется statsmodels. Используйте три каталога ранних циклов и отдельный каталог 2011–2014. Следующий шаг восстанавливает только локальный кэш идентификаторов для проверки непересечения, без переобучения или изменения замороженной модели. `LOCAL12` должен находиться вне репозитория.

```bash
export OLD_CDC=/path/to/inputs1999_2004
export MID_CDC=/path/to/inputs2005_2008
export LAST_CDC=/path/to/inputs2009_2010
export NEW_CDC=/path/to/inputs2011_2014
export LOCAL12=/path/outside/repository/replay12_inputs
python3 - <<'PY'
import os, json, shutil
from pathlib import Path
import numpy as np
from research.mortality.transport_no_crp12 import cohort, ROUTINE, sha
from research.mortality.replay_transport12 import MODEL_PATH, MODEL_SHA
roots = {y: Path(os.environ['OLD_CDC' if y < 2005 else 'MID_CDC' if y < 2009 else 'LAST_CDC'])
         for y in (1999, 2001, 2003, 2005, 2007, 2009)}
d, _, used, _ = cohort(roots)
expected = json.loads((MODEL_PATH.parent/'expected_development_sources.json').read_text())
if len(expected) != 87 or any(used.get(r['file'], {}).get('sha256') != r['sha256'] for r in expected):
    raise ValueError('Historical source fingerprints differ')
d = d[np.isfinite(d[ROUTINE]).all(axis=1)]
if len(d) != 15280 or d.SEQN.duplicated().any() or sha(MODEL_PATH) != MODEL_SHA:
    raise ValueError('Development population or model differs')
out = Path(os.environ['LOCAL12']); out.mkdir(parents=True, exist_ok=False)
shutil.copyfile(MODEL_PATH, out/'model_bundle12.json')
shutil.copyfile(MODEL_PATH.parent/'model_sensitivity12.json', out/'model_sensitivity12.json')
np.save(out/'LOCAL_ONLY_DEVELOPMENT_IDS.npy', d.SEQN.to_numpy())
(out/'development_sources.json').write_text(json.dumps(list(used.values()), indent=2))
print('Prepared frozen-model replay; participant IDs remain local, not publication data')
PY
python3 -m research.mortality.transport_no_crp12 evaluate \
  --source "$NEW_CDC" --trained "$LOCAL12" --out replay12_results \
  --expected-model-sha 40069d1a98a29308bee26e5dac895388f4c9008b983bc834ad8ac0a57fbcd094
python3 -m research.mortality.sensitivity_assay12 evaluate \
  --source "$NEW_CDC" --trained "$LOCAL12" --primary "$LOCAL12" \
  --out replay12_sensitivity \
  --expected-sha 7c7dc918f042ff44f141a9ad3b49e6573db7f9732c5eb9d93e29b09b0686fa24
```

Не перезаписывайте существующие результаты. Полное повторное обучение доступно через команду `train`, но новый timestamp закономерно изменяет файловый SHA. Это не основание менять исходную запись фиксации; для точного повторения опубликованных прогнозов используется именно сохранённая модель, как выше.

## Состав поставки и проверки

В git опубликованы научный код, обе модели, протоколы, компактные результаты, тесты и отчёт. Сопровождающий архив содержит полные агрегированные метрики, интервалы, календарную кросс-проверку, половозрастные подгруппы, аудит источников и Excel. Индивидуальные строки, идентификаторы, прогнозы участников и исходные XPORT/mortality-файлы не распространяются.

203 локальных теста выполнены без ошибок/пропусков: 64 новых, 92 Workbench и 47 реестра доказательств. Это выбранный набор, не полный повтор всех исторических научных тестов. Для четырёх тестов новых индивидуальных данных нужен `TRANSPORT12_SOURCE`; без него они пропускаются, а не считаются выполненными. CI запускает автономные проверки без данных участников; его фактический результат нужно проверять для конкретного commit.

Все интервалы условны относительно фиксированных моделей и выбранной регуляризации. Для AUC с неопределёнными повторениями отдельно сохраняется их число и ограничение интерпретации. Автоматическая публикация или развёртывание клинического сервиса не выполняется.
