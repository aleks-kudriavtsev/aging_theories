# NHANES III public cohort feasibility — 2026-10-29

Read-only metadata and history audit. No individual records downloaded, no agreements accepted, no repository files changed.

## Decision

NHANES III1988–1994 is a real publicly released long-followup route, not a controlled-access-only dead end. The public2019 mortality file is listed officially, giving roughly25–31years of calendar followup depending on baseline. An8-year event-count/ascertainment qualification is feasible, but **full frozen compact4/expanded13 prediction is not yet qualified** because the exact cardiovascular-history definition and historical assay calibration are unresolved. A feasible public data inventory is not a completed external validation.

Suggested next real computation: join first-exam cystatinC/B2M/BTP to original adult demographics/design and public mortality, report joint counts, age40–75 eligibility, assay-numeric availability, relevant survey-weight coverage, events by8years and survivor followup. No clinical prediction until gate solved. Do not require new registration or use HRS/UKB records for this route.

## Exact metadata and retrieval route

Official release index:
https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx

Verified metadata URLs:
- Adult fixed-width SAS layout: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/adult.sas
- Adult documentation/codes: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/ADULT-acc.pdf
- Main lab SAS labels: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/lab.sas
- Original cystatinC codebook: https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.pdf
- First-exam B2M/BTP: https://wwwn.cdc.gov/nchs/data/nhanes3/36a/SSNH3BTP.htm
- Repeat convenience-sample warning: https://wwwn.cdc.gov/nchs/data/nhanes3/37a/SSNH3CYS.htm
- Mortality directory: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
- Read-spec metadata listed by directory: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/R_ReadInProgramAllSurveys.R

Exact mortality URL verified from directory listing (not fetched):
https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_III_MORT_2019_PUBLIC.dat
Published bytes1663910.

Data basenames specified by official release/codebooks:
- 1a/ADULT.DAT; expected HTTPS endpoint https://wwwn.cdc.gov/nchs/data/nhanes3/1a/adult.dat
- 27a/SSCYSTAT.XPT; expected https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.xpt
- 36a/SSNH3BTP.xpt; expected https://wwwn.cdc.gov/nchs/data/nhanes3/36a/SSNH3BTP.xpt

Raw data hrefs were NOT independently resolved in this pass: web exposes numbered links but fetching a numbered data link would retrieve individual records. The local public-metadata request returned network403. Before retrieval, an authorized runner should fetch the release-index HTML and resolve actual hrefs matching /1a/+adult.dat, /27a/+.xpt and /36a/+.xpt. Avoid assuming case or promoting these expected URLs to already-verified downloads.

Machine-readable expected/verified source map and fixed-width positions: /tmp/public27_sources.json.

## Minimal adult layout

Adult SAS record content width3346,3348includingCRLF. One-based inclusive positions:
SEQN1–5; HSSEX15; HSAGEIR18–19; HSAGEU20; HSAITMOR21–24;
SDPPHASE42; SDPPSU6 43; SDPSTRA6 44–45; WTPFEX6 61–69;
WTPFHX6 70–78; MXPAXTMR1237–1240.

Adult PDF confirms HSAGEU=2(years). HSAGEIR is interview/screener age, topcoded90+; official notes recommend HSAGEIR+HSAGEU for most analyses. MXPAXTMR is MEC-exam age in months and can differ after a birthday. Either is explicit scientific choice; do not call HSAGEIR the exact age at specimen collection. For qualification age40–75, choose before reading data and retain origin in output.

Use WTPFEX6 for original MEC domain, SDPSTRA6/SDPPSU6 for combined1988–94 design. For the first-exam cystatin-selected protein domain use **WTCYPEX6**, supplied in both protein files. Check duplicate weight agreement bySEQN; do not replace it by ordinary MEC weights. WTCYPHX6 corresponds to MEC+home, a different domain.

## New markers and assay gates

1. Original first-exam cystatinC: **SSCYSTAT,mg/L**, measured with DadeBehringN-Latex/PENIA on stored baseline serum. Sampling includes all60+ with specimen,25%random12–59 and high-creatinine supplement. Codebook gives7670nonmissing/7721records, not our joint adult count. Locator fields explicitly say access/use constraintsNone. Original2006assay should not be described as ERM-DA471-standardized without applying a separately documented calibration.
2. First-exam **SSB2M,mg/L** and **SSBTP,mg/L**, SiemensN-Latex assays on stored baseline sera measured2009. Official fields report7705/7700nonmissing among7807records. These are marginal counts. Weight adjustment inherits cystatin selection and demographic/creatinine-stratum adjustment.
3. Do not use **37a/SSNH3CYS** as baseline: SSCYSCR/SSB2MR/SSBTPR are repeated second-exam values from876convenience records, with no population sample weights. Collapsing repeat/first exam would misalign baseline and sampling.
4. First-exam protein pages lack confirmed per-person LOD/LOQ/PASS flags. Positive numeric means reported-measured only, not proven QC-pass or aboveLOQ. SSNH3BTP prose gives a generic assay range inmg/dL while exact variable labels and frequencies aremg/L: preserve released field units and flag inconsistency, do not multiply by10 or screen using guessed limits.
5. Cystatin calibration drift across2006/2009 is explicitly documented. Published primary calibration study: https://pmc.ncbi.nlm.nih.gov/articles/PMC3703771/ . For qualification leave published measurements labelled original; do not silently harmonize.

## Frozen panel feasibility

Verified main-lab candidates (field/unit): CEP serumcreatinine mg/dL; AMP albumin g/dL; GHP HbA1c%; UBP urinealbumin ug/mL; URP urinecreatinine mg/dL; TCP totalcholesterol mg/dL; HDP HDLmg/dL; GGPSI GGTU/L; APPSI ALPU/L; TBP bilirubinmg/dL. Thus GGT/HbA1c/UACR are present, not absence blockers. UACR candidate100*UBP/URP mg/g requires release missing-code/QC decoding. Additional laboratory markers CRPmg/dL, FBPfibrinogenmg/dL, FRPferritinng/mL exist. CRP0.21 is released below-detection code and is not a measured high-sensitivity CRP concentration. FBP valid marginal9350, not joint availability.

Clinical exact-definition blocker: adultSAS confirms doctorCHF(HAC1C),doctorstroke(HAC1D),doctorMI(HAF10). It does not establish equivalents of modern doctorCHD(MCQ160C) and doctorangina(MCQ160D). HAF1–9 are Rose chest-pain symptoms. They cannot substitute for diagnosedCHD/angina and cannot justify setting unknown disease to0. Full frozen cvd_history is OR across5conditions; all3confirmed negatives do not establish all5negative. Status should remain unresolved until equivalent released fields are identified or a new harmonized core is explicitly developed and frozen elsewhere.

BPmed candidate HAE5A; diabetes HAD1/HAD4(gestational distinction); skin/othercancer HAC1N/HAC1O require published skip/code harmonization. No numeric risk coefficients should run on assumed normal missing states.

Historical laboratory transport is another gate: original NHANESIII creatinine requires standardization per Selvin2007, unlike continuous2001–04. Primary source https://pubmed.ncbi.nlm.nih.gov/18037092/ ; fulltext final equation reportedly −0.184+0.960*CEP, but the numerical equation was not independently recovered from official CDC metadata in this pass and should not be activated by this memo. Official lab usage notes also warn stored-serum freeze/thaw may degrade some chemistry analytes: https://wwwn.cdc.gov/nchs/nhanes/nhanes3/Lab1UsageNotes.aspx .

## Prior-use history and scope

See /tmp/public27_history.md for all ref/SHA snapshots and searches. No prior NHANESIII individual loader/results identified for the project's empirical frozen1999+ models. However multicause20 already considered III and tried archived SAS/readme metadata; completed data analysis used NHEFS. aging_biomarkers embeds author-published PhenoAge trained onIII and literary B2M/RDW coefficients derived fromIII. Therefore do not callIII wholly unknown to programme, or an independent validation of original PhenoAge/literature coefficients. A new person-level qualification/test of frozen1999+ empirical packages has a narrower defensible novelty claim.

## Remaining limits

No source file bytes, source hashes, joint participant counts or mortality events were computed here. New raw XPT/DAT retrieval is pending parent-approved pipeline, not claimed successful. The exact frozen clinical mapping is scientifically blocked, independent of legal public access. No account/agreement/action by user is needed merely to retrieve this public qualification subset.
