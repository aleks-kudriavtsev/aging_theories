# Registered national causes, iteration10

This module harmonizes national registered underlying causes from WHO Mortality Database, not modelled WHO GHE or person-level NHANES. Scientific report: `../ITERATION_10_RU.md`.

The common observed year in the retrieved standard files is2019. Russia reports list101; Germany and the USA list104. Do not describe this as the globally newest national release. The additional WHO archive returned HTTP403 and is not included.

The65-group partition exactly reconciles deaths by country, sex and age. Thirty-five groups are operationally compatible with the target natural noninfectious definition. `top20_counts.csv` contains twenty groups per country within that identified subset, not the final ranking of all target deaths. Definitions/labels are in the code's POLICY and generated `policy.json`. Bounds on ranks are conditional algebraic bounds at fixed grouping, not confidence intervals.

Reproduce with seven official archives and the original manifest:

```bash
python -m research.mortality.registered_causes10 --source WHO_MDB_INPUTS --out NEW_RESULTS
REGISTERED10_SOURCE_DIR=WHO_MDB_INPUTS python -m unittest discover -s tests -p 'test_registered10.py' -v
```

Standard library Python3.10+. Original download links are in `.github/workflows/registered-causes10-inputs.yml`; fixed hashes are in the module. Required files must match both the manifest and frozen scientific hashes. No output overwrite is permitted.

Generated outputs include585 harmonized rows,11115 sex-age rows,60 ranked rows,15 identification bounds and26 model-head coverage contracts. Full derived CSVs, Excel with formulas and the manifest accompany the downloadable iteration10 package. Raw WHO archives are not redistributed. Obtain them from https://www.who.int/data/data-collection-tools/who-mortality-database . WHO is the original data source, not the author or endorser of these analyses. Source terms require noncommercial research use; a commercial deployment licence is not implied.

USA has no2019 age-sex population in the retrieved population file. Only the total2019 Census population is used for its overall crude rate. Missing age-specific rates remain missing. No age-standardized rates or individual-country baselines are invented.

No new individual-risk coefficients or AUC improvement are claimed. Workbench and competing09 models remain unchanged. The44 new local tests included real source files; remote CI deliberately runs only36 source-free coding tests. Previous full311 scientific tests were not all rerun; see TEST_STATUS_10.json.
