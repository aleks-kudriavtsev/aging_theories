# Iteration26: frozen panels in NHANES2017–2018, one year

2026-10-29. Base aging_theories2980746f46e884ca124ccb66cc6be8f136328fb2.
Written before first individual2017–2018 input retrieval or metric inspection.
Panel25 results are known. No further urate forms or coefficients are tried.

## Why this cohort

The checked remote branch histories contain NHANES1999–2016 analyses and
explicitly exclude2017–2018 retrieval in compact14/ANALYSIS_LOCK.md. No prior
DEMO_J/BIOPRO_J/mortality2017–2018 analysis found in those refs. This supports
new-to-the-checked-history temporal assessment, not proof of no prior access
in unavailable private archives and not independent geographical validation.
CDC public mortality is through2019;2022 linkage is restricted. Use only public
2017–2018 survey +2019 LMF; no registration/agreement required. Primary horizon
one year; do not extrapolate a four-year assessment from these data.

## Fixed scientific comparison

Age40–75, linkage-eligible with valid status/exam followup, positive finite
WTMEC2YR, complete frozen predictors. Retain early administrative censors.
Primary: binary all-cause Brier1y difference expanded13 minus compact4 on the
same complete expanded13/compact4 domain, with prespecified CDC assay bridge.
Report each Brier, AUC, observed and predicted incidence, bias and O/E; paired
AUC difference secondary. Panels have different historical development and
number of cause components: this is comparison of two frozen model packages,
not controlled attribution of incremental value to the extra biomarkers.
No coefficient estimation, update, threshold selection or recalibration.
The urate coefficient from25 is not activated or tested here.

Expanded13 exact frozen_expanded13.json SHA87a8623b9fd5cb793ce2da5311a17a56b7d94251da5df335de3d406e437a55da,
developed1999–2006 per inherited archive provenance (parent archive absent in
this checkout). Compact4 exact model_compact14.json hash stored in LOCK.json,
developed1999–2010. Neither may be selected/tuned using2017–2018 outcomes.
Use frozen transformations/rates and total hazard for binary all-cause risk;
no detailed cause-performance claims from coarse late-public mortality data.

## Laboratory bridge fixed from official metadata

BIOPRO_J changed to Roche Cobas6000. Primary: published backward equations
map to DxC660i (CDC says no adjustment needed DxC800 versus660i):
albumin LBXSAL=1.044*x+0.01128 (g/dL);
ALP LBXSAPSI=10**(1.001*log10(x)-0.04294) (U/L);
creatinine LBXSCR=1.051*x-0.06945 (mg/dL);
GGT LBXSGTSI=0.8042*x+2.363 (U/L).
Bilirubin LBXSTB adjustment not recommended. No adjustment fitted here.
Nonpositive source/corrected values become missing, never arbitrary floors.
TCHOL/HDL use their own reference-method fields, not BIOPRO cholesterol.
CBC reported unchanged2015–2016→2017–2018. Inspect other required official
codebooks before processing; unresolved variable/method identity stops route.
Secondary sensitivity: unbridged released values on EXACTLY primary common
participants; retain both results, no choice by performance. Bridge uncertainty
and historical assay drift not represented by the performance interval.

## Censoring, uncertainty and quality gates

Use existing audited reverse-KM IPCW ipcw14 for all analyses, whether or not
anyone censored before1y. Deaths at horizon are events, censor at horizon known;
zero months assigned0.5month as before. Do not drop early censors or call them
survivors. IPCW assumes adequate marginal independent censoring; fail if no
horizon support. Cross-check weighted IPCW incidence with independent weighted
Aalen–Johansen. Public followup times may have synthetic perturbations.
1000paired rescaled-PSU bootstrap seed20261030, draw in all eligible MEC domain,
restrict to fixed common complete profiles, refit censoring distribution in
replicates. Models remain fixed. Count undefined replicates and do not label
percentiles among defined-only draws nominal95% coverage if failures occur.
No baseline development, bridge-estimation or model-selection uncertainty.

Report cohort flow, missingness, compact-only coverage (no cross-domain
performance comparison), weighted effective n, PSU/strata/event counts,
minimum survivor followup, censoring, per-feature winsorization frequency.
Verify SHA/byte lengths, official source URLs, within-cycle one-to-one IDs,
mortality with two readers and absence of2013/2015 participant overlap by
survey ID bounds. Preserve all result-bearing imports via source-tree file
hashes. Publish aggregate provenance and source manifest in bounded stdout
as well as artifact; no participant IDs, rows or predictions. Preserve old
results and draft PR17/118. New work in separate branch/draft PR, no merge.

Sources read before data: https://www.cdc.gov/nchs/linked-data/mortality-files/index.html
https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm
https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BPX_J.htm
https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/CBC_J.htm
