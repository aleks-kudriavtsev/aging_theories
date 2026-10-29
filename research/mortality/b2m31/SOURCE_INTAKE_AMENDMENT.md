# Source-intake amendment before any comparative result

Initial protocol6774fbe and scientific implementationdb36608 were followed by
CI run40907895237, empirical job123567576780. Intake retrieved50 files, but
the two www.cdc.gov mortality-documentation PDFs returned HTTP403:

- https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf
- https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf

The intake failed; the comparison step was SKIPPED. No participant parser,
model, tuning or new Brier/AUC/calibration comparison ran in that job. The50
downloads were not yet passed through historical-content checks. We do not
claim that their data hashes matched merely because retrieval succeeded.
No denied PDF request is retried and no alternate-host/ZIP route is used.

The endpoint definition/layout was already reviewed from public NCHS metadata
and is implemented in two pinned independent readers. The three mortality
ASCII sources are exact previously used2019 files whose historical SHA must
match before parsing. Therefore a fresh network copy of these two background
PDFs is not necessary to establish a new endpoint or a new release here.
The retained `validation08/reader_reference_audit.json` binds the official CDC
R specification (SHA04c514fd6b798ceb7eeed1afea61f42c6cdd270734c1483e1803027ff9066ac3)
and its historical exact-reader check of31126rows/249008cells for these sources.
That audit and cohort27 semantic review are additionally pinned in LOCK.
The former runtime-PDF requirement is narrowed ONLY for these two identified
URLs, with an explicit denied/not-retried provenance record in every manifest.
All other50 required inputs, selected fields, published assay marginals,
official URLs, data fingerprints, file-byte checks, full design inventory and
two-reader exact agreement remain mandatory. Unknown/missing data or any other
failed metadata source remains blocking. This is not a relaxation of readers
or permission to process an unverified mortality release.

This amendment is recorded before any comparison. It changes no population,
weight, endpoint, horizon, predictor, transformation, fold, penalty, contrast,
sensitivity or uncertainty rule. The same initial plan remains pinned; LOCK
additionally binds this amendment. Reader disagreement errors are also wrapped
to suppress accidental row-value diagnostics while still stopping the run.

Independent review and the original failed run are preserved with final results.
The two PDFs were not newly downloaded or independently rehashed in this
iteration; do not claim otherwise. Current source availability, frozen historical
layout and valid numerical comparison are separate evidence questions.
