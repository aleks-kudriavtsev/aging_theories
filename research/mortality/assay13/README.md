# Assay robustness 13

Post-validation sensitivity on already examined NHANES2011-2014 profiles; NOT new independent validation or measured assay bias. No model coefficients or default Workbench behaviour changed. All45 engineering scenarios retained for both frozen models.

- [Quantitative report](REPORT_RU.md)
- [Fixed plan](ANALYSIS_PLAN.md), commit95458703478ff686ad29b0ae20bcd4430e682309
- [Paired laboratory study protocol](PAIRED_STUDY_RU.md)
- [Selected test status](TEST_STATUS.json)
- [Historical baseline provenance](BASELINE_REFERENCE.json)
- [Metadata template](metadata_template.json) and [CSV header](pairs_template.csv)

The full companion archive includes180 scenario rows,900 cause-shift rows,64 subgroup rows,48 interval rows, source fingerprints and the detailed ITERATION_13_RU.md. Participant records are not republished. Scientific output is reproduced using the scripts, not reconstructed from rounded report numbers.

Use the exact previous transport12-evaluation1114 artifact from Actions run34456564137. Its retention is limited; do not assume perpetual availability. Official URLs are in .github/workflows/transport12-evaluation-inputs.yml. Input manifest SHA256 must equal0fd86c5dca2b55ec9a9a869f513ac9ae72ccb4c46020fce1f22ef618b12803d3. A different retrieval is a changed source version requiring review. BASELINE_REFERENCE.json records companion-archive provenance, not a nonexistent historical Git file.

```bash
python -m research.mortality.assay_sensitivity13 --source CDC_INPUTS --out NEW_RESULTS
python3 -S -m research.mortality.method_bridge13 --demo --out synthetic_bridge.json
python3 -S -m research.mortality.method_bridge13 metadata.json --pairs paired.csv --out comparison.json
python3 -S -m unittest discover -s tests -p 'test_method_bridge13.py' -v
ASSAY13_SOURCE_DIR=CDC_INPUTS python -m unittest discover -s tests -p 'test_*13.py' -v
```

Scientific calculations require previous numpy/pandas/scipy/scikit-learn dependencies. The paired-method module and46 tests are standard-library only. Twelve source-dependent tests skip without official inputs; do not describe that as the complete68-test execution. Metadata template blanks intentionally block use. Synthetic success, declared provenance and matching method names cannot authorize clinical deployment.
