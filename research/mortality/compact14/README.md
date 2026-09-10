# Compact-panel research14

Four reported measurements plus clinical profile: HbA1c, serum creatinine/eGFR, UACR and serum albumin. UACR requires two urine analytes; this is at least five underlying laboratory determinations, not measured cost reduction.

## Frozen sequence

- Analysis lock2389c6e before new2015 data retrieval.
- Coefficients96f0ca8 before new2015 data retrieval.
- Initial coverage check stopped before performance evaluation: four survivors had35months of follow-up.
- Ascertainment amendmentb105c62: IPCW, unchanged3year primary horizon, models and participants.
- New test common domain2922,93 deaths within3years. No participant overlap or outcome-based model adjustment.

Read REPORT_RU.md and RESULTS_SUMMARY.json together. Compact4 improves the primary criterion versus clinical; noninferiority/equivalence to six6 or full8 was not tested or established. All-cause includes infections/external causes. Late data identify three death groups; five latent heads are not each validated.

## Runtime (standard library)

```bash
python3 -S -m research.mortality.replay_compact14 --demo --out synthetic14.json
python3 -S -m research.mortality.replay_compact14 input.json --out research_result.json
python3 -S -m unittest discover -s tests -p 'test_compact14_runtime.py' -v
```

Only research US40–79years,1/3years, complete compact profile. Strict metadata and unit checks. Default Workbench and older coefficients remain unchanged. No clinical use or automatic assay correction.

## Scientific reproduction

Requires numpy,pandas,scipy,scikit-learn and the prior official development inputs. Exact training code is ../compact_panel14.py; its byte hash is stored in model_compact14.json. Training creates LOCAL_ONLY_DEVELOPMENT_IDS.npy for the disjointness audit; these identifiers must not be committed/distributed. Timestamp changes on retraining, so numerical parameters are compared separately from whole-file provenance. Use the original frozen JSON for exact published evaluation.

```bash
python -m research.mortality.compact_panel14 train --old-source OLD1999_2004 --source08 OLD2005_2008 --source09 OLD2009 --out TRAINED
# Keep generated IDs locally; evaluate with the published frozen model JSON in TRAINED.
python -m research.mortality.evaluate_compact14 --source NEW2015 --trained TRAINED --out NEW_RESULT_DIR --expected-model-sha 274753ea36728652d10b10fa7d02b6f77b37304084bb2b0b056128b734e394d4
NHANES_COMPACT14_SOURCE=NEW2015 python -m unittest discover -s tests -p 'test_compact14_science.py' -v
```

Do not call compact_panel14.py evaluate for the amended analysis: it preserves the original stop-on-early-censoring protocol. The separate amended module is required. New result directories are mandatory; scientific results are not overwritten.

The accompanying archive includes full aggregate metrics, intervals, source fingerprints, training CV, data-flow and numerical audits, extended ITERATION_14_RU.md and XLSX. Git contains the code, immutable protocols, coefficients and concise results. The archive does not include individual source data, local participant-ID arrays or private predictions. Source retrieval workflows have no schedule and no repository-write permission. Expired input artifacts require a new declared retrieval, not a silently substituted dataset.
