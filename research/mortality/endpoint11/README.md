# Endpoint feasibility 11

Retrospective analysis of previously studied NHANES 1999–2010 cycles; no new model or independent performance claim. General examination sample, age40–79:17,629 eligible and15,273 complete profiles. The old protein-subsample eligibility is not reused.

Read `../ITERATION_11_RU.md` for event counts, exact scope and limitations. `DATA_ACCESS_PACKAGE_RU.md` defines the detailed data request; no application has been sent. `results.json` is the compact aggregate result. Full event/interval/missingness tables, the87-file source manifest and workbook are in the accompanying mortality_endpoints11 archive and are reproducible below. No participant-level rows are redistributed.

```bash
python -m research.mortality.audit_endpoint_cohorts11 --old-source OLD --source08 INPUTS08 --source09 INPUTS09 --out NEW_RESULTS11
ENDPOINT11_RESULT_DIR=NEW_RESULTS11 python -m unittest discover -s tests -p 'test_endpoint11.py' -v
```

OLD and INPUTS08 are the official inputs of earlier iterations; INPUTS09 is retrieved by endpoint11-source-audit.yml. All require original manifest.json. Fresh downloads must be checked against the archived scientific file hashes; changed source data are not silently treated as the same version.

The offline contract has45 tests; six further tests use real aggregate results. Their passing is not clinical validation. Old model coefficients and Workbench are unchanged. Upcoming2011–2014 outcomes have not been inspected; CRP is not listed in their laboratory catalogues, so a reduced model must be fitted/frozen before validation.2015–2018 has only3 public cause labels, not negative stroke/renal events.

Public-source dictionary download workflow failed in this session; reference text was read online, no PDF visual inspection or PDF checksum claimed. Data XPORT/DAT retrieval succeeded. Historic outcomes cannot be reconstructed by matching public and restricted linkage releases.
