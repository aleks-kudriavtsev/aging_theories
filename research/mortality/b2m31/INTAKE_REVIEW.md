# Independent review: pre-comparison source-intake amendment31

Reviewed **2026-10-29**, against scientific commit `db366080` and failed empirical run **40907895237**, job **123567576780**. Read-only review of retrieval/cohort/endpoint code, existing primary-source reviews and the saved failed-job log. No network requests, PDF retries, alternative hosts, data-record reads, model fits or outcome comparisons were performed.

**Decision: scientifically sufficient to proceed under the narrowly specified, committed amendment. The two runtime documentation403 responses alone are not a substantive blocker to this already source-defined analysis.** They prevent claiming a fresh byte-level retrieval/verification of those PDFs. They do not erase the previously read primary definitions or change the identity of the three frozen2019 mortality inputs. This conclusion is conditional on preserving the full data/schema/reader/endpoint safeguards below; it is not approval to proceed with an unknown endpoint definition or changed data.

## What actually failed

The saved job log records `retrieved:50` and exactly two failed requests:

| Runtime target | Official URL | Observed result |
|---|---|---|
| `mortality_dictionary.pdf` | https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf | HTTP403Forbidden; **denied**; no successful bytes or SHA claimed |
| `mortality_description.pdf` | https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf | HTTP403Forbidden; **denied**; no successful bytes or SHA claimed |

The exception was raised in `retrieve()` before the workflow could run `b2m31.analyze`. The log contains no `B2M31_POINT_COMPUTED` or model-comparison output. Retrieval success for50 files is not yet proof that their scientific hashes/columns passed cohort verification; those checks must run unchanged in the fresh amended execution. The initial failures are known not to have succeeded, not pending or uncertain. Neither PDF request should be retried or routed through an alternate host, account, browser, mirror, proxy or archive as part of this amendment.

## Independent retained basis for endpoint interpretation

Before this failure, the same official2019 dictionary and description were read through the public web source and summarized in the already committed/pinned `b2m31/MEASUREMENT_REVIEW.md` and `SOURCE_WEIGHTS_REVIEW.md`. They specify:

- NHANES linkage key `SEQN`; adult public linkage eligibility `ELIGSTAT=1`.
- `MORTSTAT=0/1` means assumed alive/deceased; missing or ineligible status is not alive.
- `PERMTH_EXM` is months from MEC/exam to death or31Dec2019; `PERMTH_INT` has a different origin.
- The binary outcome is death through96 exam-origin months, with zero-month deaths retained; later deaths are event-free at96months. Eligible alive observations shorter than96months must stop the current binary analysis.
- The released follow-up/cause fields may be perturbed for disclosure protection; vital status is not perturbed. New2022 linkage is restricted and is not substituted here.

The earlier `cohort27/INDEPENDENT_SOURCE_REVIEW.md` independently records that the official description uses the same layout across surveys and confirms exam-origin follow-up. More specifically, `audit_reader_reference08.py` and **`validation08/reader_reference_audit.json`** retain a source-specification check against the official CDC R read program, SHA256 **`04c514fd6b798ceb7eeed1afea61f42c6cdd270734c1483e1803027ff9066ac3`**, NHANES section lines175–187. It was parsed as text, not executed as R.

Both present readers use the same previously reviewed1-based endpoints: SEQN1–6, ELIGSTAT15, MORTSTAT16, UCOD_LEADING17–19, DIABETES20, HYPERTEN21, PERMTH_INT43–45, PERMTH_EXM46–48. `nhanes_benchmark.read_mortality` uses pandasFWF and checks record widths/status semantics; `validate_temporal08.independent_mortality_reader` parses tokens manually. `audit_mortality_reader` requires exact equality of all eight fields, not just death counts. No reader change is needed or recommended.

The retained historical reference audit already checked the exact three proposed source hashes:

| File | SHA256 | Historical rows / cells checked |
|---|---|---|
|1999_mortality.dat|`562bd367107add0b2b40fed36a062d8b194fca31e11aab897a0efad4bebc8b26`|9965 /79720|
|2001_mortality.dat|`2788c8c8a02995ea5686882b5fe65b4b9cd415016031525f8a3c1a3646c47f1b`|11039 /88312|
|2003_mortality.dat|`639f6b3292a55b112fbfc842954af3d8d32bcce287cd8ae49319e5aca1c033c0`|10122 /80976|

These hashes exactly match the already locked `nhanes04/SOURCE_DATA_SHA256.txt`. I checked that equality and the historical sums **31126 rows /249008 cells** from repository aggregates only. These are **historical validation evidence**, not a claim that the stopped31 run completed current reader audits. The fresh amended run must reproduce the reader checks on its actually retrieved source bytes before any comparative fit.

## Required narrow amendment

1. Commit `SOURCE_INTAKE_AMENDMENT.md` and its LOCK hash **before any new comparative model results**. Preserve original plan commit6774fbe and identify the amendment commit separately; do not describe the runtime intake policy as unchanged preregistration.
2. The amendment must state the exact two URLs, first failed run/job, HTTP403 outcomes and no-retry rule. A fresh run must not claim it independently observed403 again: its status is a carried-forward documented denial, deliberately not requested in that run.
3. Keep the two PDFs in a separate, explicit unavailable-documentation section with source identity and historical denial provenance. Do not put fake content hashes, error-page hashes, zero-byte placeholders or `retrieved` status on them. Do not silently omit them from provenance.
4. Restrict the required retrieval set to the already specified **50 actual inputs/documentation files**:25 historical data files,22 component codebooks and index/weighting/variance HTML. Only these exact two runtime PDFs become non-required. A general rule tolerating all failed metadata, allPDFs, unknown files, changed source URLs or changed scientific data is not acceptable.
5. Preserve the pre-parse exact file-set/URL/size/SHA checks, all25 historical data fingerprints, the corrected SSCARD release schema/marginal audits, XPT zero audits and full DEMO/mortality key universe. The three mortality SHA checks and independent exact8-field reader agreement remain mandatory. No stale local raw fallback or replacement dataset is introduced.
6. Preserve every feature, creatinine calibration, assay/LOD policy, source cycle, weight formula, age/domain restriction, stratum/PSU inventory, endpoint/time origin/horizon, training split, tuning rule, bootstrap seed/count, primary/secondary contrasts and sensitivity policy. No semantic change is justified by these documentation delivery failures.
7. Include negative regression coverage showing that an unexpected failed/missing input or either missing mortality source still stops, and that the two known documentation entries cannot be converted into successful data. Existing reader/schema tests remain strict.
8. Record the retained scientific evidence explicitly. At minimum the already pinned measurement/source reviews and old data-fingerprint file remain pinned. Adding the historical reader-reference aggregate audit as another LOCK metadata item would strengthen inspectability without changing old results.

Useful existing evidence SHA256 values, independently computed locally:

- `validation08/reader_reference_audit.json`: `d3ba484633bc423b4ec55f692db69abff6241713f3780ccedda982d58647fb66`.
- `cohort27/INDEPENDENT_SOURCE_REVIEW.md`: `001f4c796184d83a6d5114dd6934bb11ea5f6460173c05c19b3ac53b259311f5`.
- `b2m31/MEASUREMENT_REVIEW.md`: `4ccb939bd9a5b4c2e91bfccc462cf65c737f3bd69e00759ab2a7ce270f65a874`.
- `b2m31/SOURCE_WEIGHTS_REVIEW.md`: `e645e9a47ebc3cd08642fc627bf019f088d6deaccc42414582f283d5402c870d`.

## Precise limits and stop conditions

No current runtime PDF bytes/hash were verified, so a statement such as “all52 sources successfully retrieved/byte verified” would be false. A defensible statement is “50 required inputs verified, with2 previously reviewed official reference PDFs denied403 during the first attempt and not retried; endpoint definitions and reader specification retained from the cited prior review.”

Reader agreement alone cannot establish medical meaning; its sufficiency here depends on the retained **prior primary-source semantic review plus exact historical source identity**. Conversely, a reviewer should not require a repeated delivery of identical explanatory material as though it were a new participant measurement when the actual frozen data/specification are unchanged.

If any of the25 historical data hashes differ, an input column/record layout is unresolved, the two current readers disagree, eligibility/status/time fail, or early censoring appears, the calculation must still stop. The amendment does not authorize choosing a new file, accepting restricted2022 follow-up, changing horizon, assigning unknown status to0, relaxing parsing or bypassing any HTTP403. Public-use linkage perturbation, assay uncertainty and reused-cohort interpretation remain material limitations after successful intake.
