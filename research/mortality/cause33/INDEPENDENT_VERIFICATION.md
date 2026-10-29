# Независимый аудит cause33A

**373/373 PASS; общий исход: PASS.**

Scientific commit `4ff18aa4115fa0b0d93a0addf9f5c1e0c86a8a83`; protocol `82eacda12eced3c0105019fee09236b647821f45`; CIrun40924184727,job123625177367.
Actual PR merge checkout `53002102a1af5a5c23d7df127dd137dbb0f1856b`; synthetic Python3.11/3.13: 28 tests each, completed success.

Common n=6044,allcause8y=586,unknown early=0.
Weighted burden rank: 002, 001, 003, 004, 005, 007, 008, 009, 006. Priority top5: 002, 001, 003, 004, 005. Proceed: 002, 001.

| Code | Events8 | Event PSU | Weighted events8 | Support failures | Priority | Proceed |
|---|---:|---:|---:|---:|---|---|
| 001 | 148 | 67 | 1882199.21474333 | 0 | True | True |
| 002 | 176 | 75 | 2383527.11994933 | 0 | True | True |
| 003 | 33 | 25 | 429776.74088567 | 2 | True | False |
| 004 | 25 | 21 | 422773.70826600 | 7 | True | False |
| 005 | 25 | 20 | 255918.78017700 | 7 | True | False |
| 006 | 4 | 4 | 66604.68824333 | 42 | False | False |
| 007 | 23 | 19 | 203250.32574133 | 9 | False | False |
| 008 | 13 | 13 | 126712.78814900 | 30 | False | False |
| 009 | 10 | 8 | 100497.45151733 | 33 | False | False |
| 010 | 129 | 60 | 1578766.93614467 | 0 | False | False |

Проверены83supportcomparisons/cause,40foldcells/cause,all10+unknown accounting,percycleweights/ratios,independent hash folds,weightedranking/ties/no backfill,old95hashes иsource50/data25. Нет marker-benefit или clinical/external-validation вывода.

Ограничения:

- This is a source/count/ranking/support gate only. No models, comparative prediction metrics, model bootstrap, Fine–Gray or joint integration are fitted or inferred.
- Source hashes, parsed outputs, accounting, ratios, fold maps, ranks and all support failures are independently reconciled from aggregates. No raw records or source byte downloads were read, and no second person-level cause reanalysis was performed.
- Top5 by weighted burden is fixed before support; failed priorities are not replaced. Count-based selection is explicit and does not justify confirmatory inference or exclusive marker-cause assignments.
- Unknown early causes block future comparisons and are never negative labels. Known competing deaths are noncases only for the fixed-horizon cumulative-incidence binary target, not censored or alive.
- Complete-case/LOD/linkage selection, reported-value assay/storage/calibration uncertainty and prior cohort use remain. Weights and pragmatic support thresholds do not certify power, EPV, clinical validity or national/global transport.
- Both historical mortality PDF denials remain denied_not_retried without new hashes or retries. Completed GitHub CI logs may be independently read; no raw/PDF/ZIP GET is used by this audit.
