# Iteration27: longer-followup source/endpoint/assay qualification

2026-10-29. Base798610af4f2cebd83437c27fc93e03f3fb3f67a7. Written before
first individual NHANESIII input retrieval or same-person count calculation.
This is empirical source qualification, not a mortality-model validation.
Prior26 results are known;2017–2018 will not be used for new model selection.

## Public route and scope

Use official public CDC NHANESIII1988–1994 FIRST examination data: adult
questionnaire, original cystatinC release27a, B2M/BTP36a, mortality2019.
Do not use37a repeated-examination convenience sample or its different time
origin. No login, new account, restricted linkage or agreement is accepted.
All controlled-access routes remain subject to their actual conditions;
HRS person-level records never enter an AI/unattended-cloud process.

Qualification population: adult interview age40–75, HSAGEU=2(years), original
exam mortality-eligible with known death status and nonnegative exam followup,
positive finite WTCYPEX6 (Cystatin MEC subsample weight), valid survey strata/PSU.
Age is interview age, not falsely labelled exact MEC age. First-exam biological
specimens and first-exam-linked mortality are required; laboratory assay year
2006/2009 is not a later biological baseline.

Before opening person-level inputs, check metadata gate and file manifest,
protocol binding, official URLs, byte counts and SHA256. Never execute SAS
programs or inferred code. Parse only declared fixed-width numeric positions
and independently verify selected fields against pandas.read_fwf.
Join one-to-one within survey SEQN; adult file does not include every sampled
child, so counts outside adult source are reported, not treated as foreign data.
Verify two mortality readers; keep missing/zero/negative analyte values separate.

## Prespecified measurements and endpoint counts

SSCYSTAT cystatinC/P01034, SSB2M beta2-microglobulin/P61769 and SSBTP
beta-trace protein/PTGDS/P41222, all released mg/L, surplus serum from FIRST
exam. Numeric availability is finite positive released value, NOT established
above-LOQ usability. BTP documentation has inconsistent prose range units;
do not convert/clamp by that range. No new imputation, no NPX/RFU conversion,
no TNF/receptor substitution, no derived-marker proxy.
Use WTCYPEX6, not general MEC weight, for this protein-subsample domain;
verify Cys/36a shared weights agree where both observed. Report single-marker
and joint-three counts, joint deaths through8years, effective weighted n,
strata/PSU counts, followup support and early-censor count. No Brier/AUC/HR,
training, thresholds, biomarker-form searches or confidence claims in this
qualification. Count observed deaths only and do not invent outcomes for early
censors. Horizon8years lies inside frozen-model time support but does not make
those models applicable.

## Frozen-model application gate

Old compact4/expanded13 require exact cvd_history semantics (physician CHF,
CHD, angina, MI, stroke). NHANESIII has confirmed CHF/stroke/MI fields but the
same physician CHD/angina equivalence is unconfirmed; Rose symptoms cannot be
substituted. Historical creatinine/Cys/other method transport also requires
adjudication. Missing predictors cannot be set to normal/zero or omitted from
a frozen model. Thus frozen-model evaluation remains blocked independently
of successful protein/endpoint import. Distinguish public data availability,
correct import, analytical usability and model compatibility.
A new NHANESIII-compatible development/validation core requires a separate
prospectively recorded design; qualification counts do not imply an untouched
future test population or externally validated coefficients.

## Deliverables

Aggregate qualification counts and blockers; no participant records/IDs or
predictions. Preserve source/plan/code hashes and independent read audits.
Synthetic fixtures verify software only. Never report synthetic AUC as external
validation. Preserve previous results/models/PRs; new branch/draftPR, no merge.

Sources: https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx
https://wwwn.cdc.gov/Nchs/data/nhanes3/36a/SSNH3BTP.htm
https://wwwn.cdc.gov/Nchs/data/nhanes3/37a/SSNH3CYS.htm
https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.pdf
https://www.cdc.gov/nchs/linked-data/mortality-files/index.html
