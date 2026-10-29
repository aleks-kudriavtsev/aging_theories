# Independent cohort-history / measurement audit for iteration26

2026-10-29. Read-only audit of repository content and public CDC documentation;
no individual data downloaded and no empirical calculations repeated.

## Scope and limits

Inspected aging_theories 21 remote refs plus 3 local heads, 92 reachable commits
(count at audit time); aging_biomarkers 2 local heads,425 reachable commits.
The latter checkout has no fetched remote branches: this audit does not certify
all private branches or untracked external archives. Therefore say 'not previously
analysed in inspected repository history', not 'never used anywhere'.

Searches: git grep over each remote ref, and git log --all -G
'DEMO_J|BIOPRO_J|2017_2018|2017.?2018' -- research .github. The ONLY theory-history
commit matching2017/2018 is2389c6e, the compact14 decision explicitly excluding
2017–2018 retrieval. No2017participant retrieval/model/result code found.

## Used populations (not candidates for a new untouched test)

- NHANES1999–2004: original benchmark development/recalibration;
  research/mortality/validation08/README.md:1999–2002 n4549 development;
  2003–2004 n2381 recalibration.
- NHANES2005–2008: validation08 n6049,958deaths10yr, later development in
  competing09. See validation08/ABSTRACT_IYGF2026_RU.md and
  competing09/ABSTRACT_IYGF2026_RU.md (2005–2008 development n5287).
- NHANES2009–2010: competing09 evaluation n3182;
  commit99e8c15 lock; commit1bd73c3 completed analysis.
- NHANES1999–2010: transport12/compact14 development n15280;
  ITERATION_12_RU.md and compact14/ANALYSIS_LOCK.md.
- NHANES2011–2012: transport12 evaluation n2666/130deaths5yr,
  subsequently reused assay13; ITERATION_12_RU.md, assay13/REPORT_RU.md.
- NHANES2013–2014: transport12 n3013/113deaths4yr; NfL22/joint23
  n1270/37deaths; panel25 wider MEC n2867/96deaths. All used.
- NHANES2015–2016: compact14 n2922/93deaths3yr. Used;
  compact14/REPORT_RU.md explicitly says future selection on this cycle is
  development rather than new independent validation.
- NHEFS/NHANESI: multicause20 already obtained/analysed;
  commits858c7e5(lock),bb6dd40(retrieval),0b53c23(implementation/results).
- flchain and PBC in aging_biomarkers already used:
  analysis/validation/validate_flchain.py and validate_pbc.py; these do not
  provide the full frozen NHANES panel and cannot serve as a new full-panel test.
- HRS24 is an access/metadata route, no person-level modelling. Its
  README_RU.md documents HRS AI restrictions and Sensitive Data Access Agreement.
  Not a currently available automatic extension. No new agreement was accepted.

## Important qualification: literary awareness of2017–2018 exists

aging_biomarkers HEAD:data/diagnostic_panel/reports/conference-review/evidence.json
(lines around19365 and32958), final_panel.csv and final_panel.md include a
beta-carotene publication labelled NHANES2003–2006 and2017–2018 (PMID42393351).
This is literature-derived evidence, not a retrieved participant table or
an independent model test. Therefore2017–2018 should not be called wholly
unseen in the broader research programme. It can be described as a new
participant-level temporal evaluation in the inspected model history.

## Public path confirmed

https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
currently lists NHANES_2017_2018_MORT_2019_PUBLIC.dat (449623bytes), and2015–2016.
Only read the directory index, not the participant file. With follow-up ending2019,
1yr is a reasonable metadata-based primary horizon; ascertainment/censoring
must be checked. The2022 linkage is a separate restricted/RDC route, not a
public substitute. Do not extend the horizon because of an attractive result.

## Official measurement audit before outcomes

All public codebooks below retrieved via web on2026-10-29.

- TCHOL_J: LBXTC total cholesterol mg/dL; LBDTCSI mmol/L is x.02586.
  No method/equipment/site changes in2017–2018. Use LBXTC, not BIOPRO LBXSCH;
  the biochemistry bridge for LBXSCH is not applicable to TCHOL's LBXTC.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/TCHOL_J.htm
- HDL_J: LBDHDD mg/dL, LBDHDDSI mmol/L. No method/equipment/site change.
  Codebook contains inconsistent prose ('calculated' and LBXHDD) despite
  exact released field LBDHDD/direct-HDL label; preserve LBDHDD, do not derive
  HDL from LDL/TG or replace with another field.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/HDL_J.htm
- CBC_J: WBC/platelets thousand cells/uL, RDW%, NEPCT/LYPCT%. DxH800,
  no method/equipment/site change in2017–2018; no new outcome-based correction.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/CBC_J.htm
- GHB_J: LBXGH is%. Standard HPLC;2.2% specimens retested by secondary
  ultra2 HPLC; instrumentation changed for secondary method, described as
  nonsignificant, mean relative error0.5%(range−4.0…5.2%). No equation for
  changing individual values is supplied. Retain published LBXGH and flag
  secondary-method change as transport limitation.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/GHB_J.htm
- ALB_CR_J: URXUMA ug/mL, URXUCR mg/dL; UACR=100*URXUMA/URXUCR mg/g.
  Published URDACT rounded. LLOD replacements LLOD/sqrt2 already present;
  limits0.30ug/mL and5.00mg/dL. Retain released fills/flags; do not impute anew.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/ALB_CR_J.htm
- BPX_J provides BPXSY1/2/3 in mmHg after5min seated rest, three readings;
  fourth can be taken after incomplete/interrupted attempt. Preserve frozen
  1–3 positive-reading algebra; do not quietly adopt fourth-reading averages.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BPX_J.htm
- BMX_J has BMXBMI kg/m²; preserve recorded BMI. Some physical limitations
  induce missing measures; do not replace missing BMI by alternative formula.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BMX_J.htm
- SMQ_J: SMQ020=100lifetime cigarettes(1yes/2no); SMQ040=1everyday/2some
  days/3notatall. Skip after020no means missing040 is not missing smoking
  status; existing frozen algebra handles this. WTMEC2YR when joined to labs.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/SMQ_J.htm
- BPQ_J: BPQ020hypertension1yes/2no; BPQ040Aever prescribed meds1/2;
  BPQ050Acurrently taking1/2.020no and040Ano skip050A; do not classify these
  skips as unknown treatment. Use frozen skip-aware mapping.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BPQ_J.htm
- DIQ_J: DIQ0101yes,2no,3borderline,7refused,9unknown. Core questions
  described as similar1999–2016. Preserve frozen borderline→0; expose that
  definition rather than recoding after outcome review.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DIQ_J.htm
- MCQ_J contains MCQ160b/c/d/e/f (CHF/CHD/angina/heartattack/stroke) and
  MCQ220cancer history. SAS names may appear lowercase in webpage; checked
  reader normalizes column names. Existing definition combines only these
  five CVD fields, not generic 'anyheartdisease'.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/MCQ_J.htm
- DEMO_J has WTMEC2YR, SDMVSTRA, SDMVPSU, RIAGENDR, RIDAGEYR. Design fields
  are cycle-specific; prefix stratum/cluster labels when combining cycles.
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.htm

BIOPRO_J bridging was independently surfaced before outcomes and already
verified by parent. Source reports248paired2016serum samples;
CK excludes1outlier only. For frozen predictors, backward Cobas→DxC equations
albumin1.044x+.01128, ALP10**(1.001log10x−.04294),creatinine1.051x−.06945,
GGT.8042x+2.363. Bilirubin adjustment not recommended. Equations are method
bridges, not individual exact equivalence or outcome-fitted calibration;
uncertainty of bridge parameters is not covered by fixed-model PSU bootstrap.
https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm
