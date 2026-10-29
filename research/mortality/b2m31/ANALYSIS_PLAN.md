# Iteration31: B2M beyond a compact routine core on reused SSCARD_A

Written before new comparative outcome calculations. This is exploratory
research on NHANES1999–2004, previously used by this programme. SSCARD_A was
already imported, cystatin C fitted in M4_cystatinC, and B2M loaded. B2M is
chosen in light of the NHANESIII28/29 results. Neither confirmatory significance,
optimal panel selection, an untouched holdout nor independent validation of
freeze30 is claimed. Frozen30/old models/results/main/UI remain unchanged.

## Sources, target and comparison domain

Biological baseline is the original NHANES1999–2004 exam, not serum assay
2018–2020 or release2021/2022. Population: screening age RIDAGEYR40–75, recorded
sex1/2, eligible linked2019 public mortality, finite nonnegative exam followup,
positive revised SSCARD weight, RIDSTATR=2 (MEC examined), and valid masked stratum/PSU from DEMO. Use
official July2022 corrected SSCARD_A (21890 source rows), joined one-to-one
by SEQN to the matching three cycles; no SSBNP/troponin eligibility restriction.
Known source fingerprints from nhanes04 must match all overlapping data files.
Metadata and actual source manifest are exported with URLs, SHA and byte counts.

Endpoint: all-cause death by96months from exam, including all causes. Deaths
at0months are events; no artificial time shift is needed for binary scoring.
Require every eligible survivor to have >=96months recorded followup before
any comparison. Missing/ineligible vital status stays unknown. If the gate fails,
stop; do not change horizon or label censored people alive. Public followup
perturbation and retrospective use remain limitations.

Use corrected SSCARD_A weights: (2/3)*WTSSCB4Y for1999–2002 and
(1/3)*WTSSCB2Y for2003–2004, the NCHS exception for pooling the first cycles.
Do not use ordinary MEC weights or WTSSCB2Y/3 for all cycles. Masked strata
and PSU come from DEMO; cluster=(SDMVSTRA,SDMVPSU), never PSU alone. Allow
documented PSU1..3 in1999–2000. Before comparisons, require that every full-DEMO masked stratum/PSU pair is
still represented in the eligible target frame; otherwise stop rather than
silently shrink the design. A singleton eligible stratum stops the bootstrap.
Weights do not correct item missingness or selection into the common domain.

Primary common domain for ALL FOUR models: finite11routine terms, positive
finite SSCYST/SSB2M, both comment flags exactly0 (measurable). Both proteins are
required for every arm to prevent changing the comparison population. Report
the flow, eligible/common n and deaths, overlap of exclusion reasons, weighted
coverage, effective n, cycle counts and source/flag/pristine counts. This domain
differs from old nhanes04 (age40–79,10years,additional cardiac assays and other
features/weights); do not compare old absolute scores as an improvement claim.

## Compact core chosen by definitions/availability, not new outcomes

Eleven fixed terms, in order: (screening age−60)/10, square of that term, male,
current and former smoking, (SBP−130)/20, (BMI−27)/5, HbA1c%,
min(eGFRcr,60)/15, max(eGFRcr−60,0)/30, log2(UACR).
Smoking: SMQ020=2 never; SMQ020=1 and SMQ0401/2 current,0403 former;
unknown/refused categories remain missing. SBP is the research mean of finite
positive BPXSY1..3 (>=1 reading), not BPXSAR/released traditional average;
BPXSY4 is not substituted. BMI=BMXBMI; HbA1c=LBXGH original released percent,
no invented cross-over regression. No diabetes history or CHD/angina proxies.

Serum creatinine:1999LAB18 LBXSCR mg/dL, SCr_standard=1.013*raw+0.147;
2001L40_B LBDSCR and2003L40_C LBXSCR have no additional CDC correction.
Verify exact current columns; no NHANESIII .960x−.184 equation. Use the fixed
2021 race-free creatinine-only CKD-EPI equation; no Cys in the comparator and
no additional raw-creatinine term. Positive URXUMA ug/mL / URXUCR mg/dL *100
gives UACR mg/g. No NHANESIII urinary floor10. Retain released routine values
and documented substitutions; missing/nonpositive measurements stay missing.
eGFR is an estimate, not measured GFR or proof of individual calibration accuracy.

## Assay policy and arms

SSCYST/SSB2M are reported serum mg/L, measured on Siemens Dimension Vista1500.
Measurable flags0 are analytical detection-range flags, not individual QC or
proof of modern reference-material traceability. Release-specific Cys IFCC/ERM
standardization is not confirmed; do not apply older SSCYPC/III calibration or
claim transport to frozen30. Primary proteins are log2(reported values), without
imputation, outcome-selected thresholds, clipping or rank transform.

Four fixed arms: routine, routine_b2m, routine_cys, routine_cys_b2m.
**Primary contrast:** routine_b2m minus routine, primary metric delta weighted
Brier at8years. AUC and calibration are accompanying diagnostics.
**Secondary prespecified control:** routine_cys_b2m minus routine_cys, to ask
whether B2M adds information beyond another renal protein already studied here.
It is secondary/exploratory, not a second independent validation. No selection
of a winning protein, arm, subgroup, horizon or transform; no multiplicity-free
confirmatory claim or optimization of the panel.

## Training, scoring and uncertainty

Use a new weighted binary ridge-logistic model, not the old cause-specific
engine or transported frozen coefficients. Five whole-masked-PSU outer folds,
hash salt b2m31:outer assigned on eligible design before common-case filtering.
Three whole-PSU inner folds within each outer training, salt b2m31:inner:{fold}.
Same participants, folds and weights in every arm. Each actual training split
learns weighted mean/population SD (constant scale1); no imputation or tail fit.
Ridge lambda1/10/100 minimizes pooled weighted inner-heldout Brier, ties1e-12
prefer larger lambda. Normalize likelihood weights to that training n; leave
intercept unpenalized; preserve binary28 fit/convergence rules.

Report OOF weighted Brier/AUC, observed/predicted risk, bias, O/E,
offset calibration intercept and intercept+slope, plus calibration bins. These
diagnostics do not alter scores. Report undefined/nonconverged counts explicitly.

500 paired rescaled-PSU bootstrap draws, seed20261102, from the full eligible
design; same draw weights and fixed common masks/folds in all four arms.
Relearn all preprocessing/tuning/outer fits per draw. Export all500 paired
primary/secondary score deltas and statuses, no participant predictions.
Percentile2.5/97.5% intervals are conditional on the fixed partition, reported
measurements, chosen forms and programme history. They do not include assay
calibration/selection-history/external-transport uncertainty. If any metric is
undefined, disclose valid/undefined draws; defined-only bounds are not a full
nominal95% interval. No partition search or interval pooling.

## One predefined limited sensitivity

Recompute all four arms on their own shared reported-values domain: complete
routine core, positive finite proteins and known flags0/1 only. This retains
published LLOD/sqrt2 fills without calling them precise quantification.
Above-ULOD flag2 remains excluded: its exact released fill/censoring rule is
not separately documented, so this sensitivity does not treat it as a value. Report domain change and within-domain
point metrics/contrasts only; no extra bootstrap or absolute cross-domain
superiority claim. Primary retains both pristine/non-pristine samples; report
SSPRIS distribution. Glycated-albumin storage warnings are not evidence that
all non-pristine B2M/Cys are analytically invalid. No additional subgroup search.

## Delivery and stop rules

Protocol/metadata/source specs committed before empirics. Reader keeps exact
XPORT zero repair from source bytes (critical for flag0), never epsilon rounding
or relaxed validation. Check source schema, unique IDs, cross-cycle joins,
official marginal release counts, units, source hashes and mortality reader.
If official corrected weights/assay definitions/data are unusable or inaccessible,
stop the dependent analysis and report the exact gate; no synthetic validation.
Only aggregate outputs/source manifests/fold coefficients/bootstrap deltas are
published. Test source/feature/domain/nesting/weight/undefined behavior; audit
results independently and run scientific CI. Separate draft PRs, no merge,
approval, deployment, accounts, agreements, restricted requests or contacts.
