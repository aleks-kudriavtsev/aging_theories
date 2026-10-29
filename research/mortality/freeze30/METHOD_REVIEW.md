# Independent scientific method review — freeze30

Reviewed 2026-10-29; scientific code commit `16e2eb8fe0abd51b351a6c3775a7180ebdd25354`, plan commit `5cbf90e`. This is a read-only review by an agent that did not author freeze30. No participant data, archive downloads or model refits were used. Review covered adapter.py, freeze.py, qualification.py, CONTRACT.json, LOCK.json, FREEZE_PLAN.md, SOURCE_POLICY.md, the qualification reports/manifests and synthetic test definitions. The actual artifact has a separate independent aggregate/numerical audit.

**Finding: no critical scientific error found that blocks this research freeze.** It constructs two future-use research artifacts from a previously assessed cohort; it does not establish new external performance. This assessment is not a PR approval, authorization of a new data source, or clinical endorsement.

## Training, scoring and numerical scope

- Both arms use the same routine29 complete-profile domain, fixed expected n3730/events563, all three positive reported proteins, WTCYPEX6 and the preserved 96-month exam-origin endpoint. The source loader preserves the historical common-domain and censoring guards. Final training separately refuses early-censored alive participants.
- Final lambda selection uses exactly three whole-PSU folds on the full known training domain, one prespecified salt, grid1/10/100, pooled weighted inner-heldout Brier and the larger-lambda tie rule1e-12. Each inner preprocessing fit sees only its training split. The final scaler and ridge model are then learned on all training records. This is a defensible finalization algorithm and is explicitly different from claiming the inner tuning criterion as validation.
- The unpenalized intercept, weighted population mean/SD, SD1 for constants and training-n normalized ridge weights agree with the preserved numerical primitives. Both arms share participants/folds. No single-protein winner, new transformation, subgroup or outer partition is selected during this freeze.
- Prediction applies the frozen transform and coefficients without refitting, recalibration, local scaling or imputation. The routine arm still requires the complete joint input domain, although protein values do not affect its numerical score. This matched-domain restriction is deliberate and must remain visible when reporting comparisons.
- The artifact is a separate eight-year weighted ridge-logistic model, not the earlier five-cause piecewise-constant engine. An 8-year probability does not imply clinical transportability.

## Measurement and time contracts

- The exact eleven routine terms plus three log2 proteins match routine29 and the written plan. No extra creatinine or diabetes-history predictor duplicates the renal/HbA1c terms.
- Input serum creatinine is explicitly original historical NHANESIII CEP. The `.960*CEP-.184` correction is applied once, only to serum; UACR uses the separate urine measurements and the correct `100*UBP/URP` factor. Modern IDMS creatinine must not be mislabelled CEP.
- Adult2021 creatinine-only eGFR uses the reviewed sex/age equation without race or cystatin. It remains an estimate from historically calibrated material; historical calibration uncertainty is not removed by an equation.
- The adapter rejects unknown mapping/unit labels, original sentinels, null/nonfinite values, boolean numerics, unsupported age/categories and URP<10. It retains UBP0.4, released BMI edits and HbA1c method4 in the primary representation. These are documented reported-value policies, not individual QC/LOQ certification.
- Only the original first-exam reported mapping is accepted. Modern standardized cystatin, NPX/RFU, repeat examinations and other cohorts need separately reviewed mappings. No approved external mapping exists. Declaring the correct mapping string cannot authenticate actual specimen provenance, and the documentation correctly says so.

## Provenance and determinism

- The plan and measurement/source policies were committed before final fitting. The prior29 result, source manifest, metadata and all80 historical scientific Python hashes are checked; frozen construction appends new files instead of altering prior science.
- The source policy permits only a dynamic catalogue-index byte change and requires exact hashes for all13 non-index files plus identical original/resolved release URLs. This is not a new assay release or relaxation of record-reader validation.
- The two-run workflow compares full artifact and freeze-manifest bytes in one pinned environment. This supports deterministic replay in that environment; it is not evidence of platform-invariant floating-point coefficients.
- Raw measurements, IDs, outcomes and individual fitted/heldout probabilities are not exported. Training replay exports only n and maximum absolute score difference. Inner Brier/objective values are training diagnostics and must never be advertised as new Brier/AUC performance.

## Access, overlap and qualification claims

- The external-mapping registry is empty, so metadata declarations cannot enable scoring. The gate has no participant/network reader, distinguishes exact visits/assays/time origin, blocks short/unknown survivor follow-up and retains HRS-specific AI constraints. Self-attestations are explicitly not legal or source authentication.
- The actual12 candidates remain blocked. The reports correctly separate an unapproved controlled-data route from a failed metadata GET or an unresolved URL. No new independent validation or participant release is claimed.
- ARIC visit4's missing documented same-visit HbA1c, MESA-Kidney's late2012–14 baseline/2020-package limitation, HRS policy/follow-up limitations, and UKB NPX/current application pause are material scientific/access gates. Distinct cohort recruitment does not itself prove zero individual overlap.
- The SSCARD source is explicitly recognized as previously loaded/analysed in the programme, with Cys already fitted. A future B2M/two-protein question on it would be reused-cohort research, not full3 independent validation; BTP must not be imputed or replaced.

## Nonblocking limits to retain

The adapter's numerical schema validation and paired manifest hash detect errors against a reviewed reference; they do not authenticate an attacker-replaced artifact/manifest pair or independently verify the truth of provenance fields. Full source-byte verification and refitting are not available to this read-only aggregate reviewer. Real measurements cannot be analytically validated from positive numeric availability alone. Future external work still requires an approved mapping, actual endpoint/censor audit and lawful access in the actual operator environment.

The protocol does not include training uncertainty/assay uncertainty in a final coefficient interval; it freezes point parameters for a future evaluation. Old routine29 intervals remain conditional internal-CV estimates and must not be relabelled uncertainty of independent external performance of this final fit.
