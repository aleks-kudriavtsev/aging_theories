# Period39 — official Seattle Heart Failure Model source qualification

2026-10-30. Bounded review of one postdiagnosis model. No account, agreement acceptance, role selection, form submission, numerical oracle query, patient data, model predictions, fit, model code, ZIP, mirror, repository edit, publication or holder contact.

**Result: the updated2013 SHFM is not yet qualified as a fully reproducible1–5year source component.** Its official update supplies a one-year survival/score relation but not the complete updated five-year baseline array, coefficient and input contract, or a versioned permissively licensed source package. This is a source qualification result, not a claim that SHFM itself lacks clinical evidence.

## Verified official sources

1. [University of Washington SHFM update page](https://depts.washington.edu/shfm/update.php), successfully read by the web tool before the later executor denial. The page separates2013 model changes from2017 JavaScript/unit-interface changes. Original five-year survival extrapolated PRAISE0–3year survival exponentially; the2013 update instead uses observed1–5year survival from a larger SHFM-D derivation population. It gives `SHFMScore=ln[ln(S1)/ln(0.9604)]` and rounded one-year pairs `0→96.0%,1→89.5%,2→74.0%,3→44.0%,4→10.8%`. Updated2–5year numerical baselines are not listed. The text also changes LVAD/ICD and treatment logic. Its statement of higher survival conflicts directionally with its example of five-year mortality increasing22%→25%; do not resolve that by assumption.
2. [Official SHFM landing page](https://depts.washington.edu/shfm/), successfully read. It describes projected survival for patients with HF, including baseline and intervention scenarios, and presents a healthcare-provider role selector. No role was selected. This page does not supply numerical age/EF eligibility bounds, an onset-date origin, or a complete parameter table.
3. [Official linked Terms of Use](https://depts.washington.edu/shfm/term.php), read by the parallel rights reviewer. See exact scope in `/tmp/period39_shfm_rights.md`. Public accessibility is not an open-source licence. The terms cover code as content and reserve reuse beyond their stated allowance; additional model-specific terms may exist. No terms were accepted.

## Ready / blocked fields

| Required field | Status | Evidence or missing item |
|---|---|---|
| Source identity and institutional chain | Ready for registry metadata | Official University of Washington domain and linked update/terms pages. |
| Endpoint family | Partially ready | HF projected survival is stated. No separate vascular/nonvascular death partition was recovered. Do not relabel this as a disease-incidence, disease-stage, or cause-specific mortality model. |
| Population | Partially ready | Existing HF, not a disease-free general population. Exact approved model-version ranges and exclusions were not recovered from the reviewed pages. No Russian validation or transport qualification follows. |
| Prediction origin | Blocked for onset bridge | A baseline assessment is described, but no verified equivalence to date of HF diagnosis was found. Time since diagnosis must not be silently set to zero. |
| Time unit and available horizon | Partially ready | Source describes survival at1–5years. The updated5year contract must be distinguished from the original extrapolation. Ten-year mortality is not qualified. |
| Updated one-year baseline relation | Ready as a published fact | Formula uses0.9604 in the denominator. Algebraically it implies `S1=0.9604^exp(SHFMScore)` for the score defined there. This is not a complete score calculator. No values were recomputed in this review. |
| Updated2–5year baselines | Blocked | Required numerical `S0(2),…,S0(5)` or equivalent integrated hazards and their version/normalization were not recovered. Do not use an exponential continuation ofS1: that would reinstate the superseded extrapolation assumption. |
| Full updated linear predictor | Blocked | Actual coefficient vector, centering/reference terms, exact transforms and validation rules are not supplied by the reviewed update/landing text. A transformation from survival back to a score is not a definition of that score from clinical inputs. |
| Units and missingness | Blocked |2017 notes SI/pounds/default options, but no complete conversion constants, allowable ranges, default values or missing-data semantics were recovered. Missing labs must not be silently filled by class defaults in a future implementation. |
| Updated device/treatment logic | Blocked | A complete executable specification, including variable ICD benefit, NYHA overrides, dose equivalences and device eligibility interactions, is missing. A constant multiplier would not reproduce a risk-dependent rule. |
| Five-year numerical oracle | Blocked | The published rounded one-year score pairs can check a score-to-survival transform only; they cannot check coefficients or2–5year survival. No calculator request was submitted. |
| Versioned code/parameter manifest | Blocked | No authenticated immutable source artifact/hash, release/commit and complete updated parameter manifest were obtained. The current live site cannot be assumed to be an unchanged2013 implementation. |
| Reuse rights | Blocked for a claimed open code port | No permissive software licence verified. Distinguish copyrightable source-code copying/translation from independently implemented mathematical facts; this review does not decide rights to an unseen complete publication or grant any permission. |
| Russian40–80,5/10year target | Blocked | Requires its own source-to-target eligibility, baseline, treatment-era and outcome validation. HF survival alone cannot identify general-population onset→mortality probabilities. |

## Mathematical implications and next gate

These are our methodological deductions, not extra parameters recovered from the website.

A Cox-style survival relationship `S(t|x)=S0(t)^exp(eta(x))` needs both the full predictor and the baseline at the requested horizon. One baseline point does not determine the remaining function: infinitely many monotone survival curves have the same one-year value and different five-year survival. Therefore the one-year facts above cannot identify an updated five-year prediction.

Combining this postdiagnosis component with an incidence component would additionally need a common disease definition, diagnosis/severity/time origin, mortality after different onset times, competing events and a calibrated Russian bridge. Multiplying a five-year incidence by a five-year postdiagnosis mortality is invalid: people develop HF at different times and have different remaining follow-up. No source-to-target coefficients or hazards should be inferred from the parent's period population-mortality reference.

A subsequent protocol could qualify a complete open mathematical specification independently of the site's source code only if it fixes the model version, all predictor/therapy rules, numerical1–5year baselines, unit mappings, domain and a reference. It must not combine original2006 coefficients with an unspecified2013 baseline and current interface logic. Until those fields close, retain the registry entry with `absolute_prediction_ready=false`; do not implement or publish a purported updated SHFM calculator from this review.

## Exact read outcome and restrictions respected

- **Successful:** web reads of official update and landing pages. Parallel reviewer read the official linked terms. These supplied source facts only; no model was executed.
- **Denied:** read-only executor `urllib.request.urlopen('https://depts.washington.edu/shfm/update.php')`, intended to inspect static HTML/script links and obtain a source-byte checksum, returned `<urlopen error Tunnel connection failed:403 Forbidden>`. This exact route was stopped; no retry, headers, mirror or alternate-transport source retrieval followed. No HTML source file was written because the read failed.
- **Pending prerequisite:** complete authenticated public non-ZIP updated baseline/coefficient/units/reference package, plus appropriate reuse qualification. This is not a pending user approval. No attempt was made to obtain it by contacting the holder or accepting a website role/agreement.
- **Uncertain operations:** none. No write or calculator-submission outcome is uncertain. Only local review notes were written.

Supporting rights note: `/tmp/period39_shfm_rights.md`. No original-code checksum is claimed.
