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

## Addendum: complete currently fetched aging_biomarkers remote-branch scan

2026-10-29, after parent fetched all remote branches: independently scanned
98 origin refs and469reachable commits. This supersedes the initial limited
aging_biomarkers scope above; it does not claim access to deleted/unreachable
commits, untracked archives or other repositories.

Programmatic git grep -n -I -E across each of98refs, restricted to analysis,
tools,data/model_evidence,.github,docs, with pattern
`DEMO_J|BIOPRO_J|2017_2018_MORT|2017.mortality|2017_mortality|2017_.*_J[.]xpt|NHANES.?2017`
returned0current-tree matches. git log --all -G for the same pattern and paths
returned0commits. Thus no person-level2017–2018 retrieval/analysis/source
provenance was identified in the fetched history.

The broader code-history pattern2017.?2018 in analysis/tools/.github returns
one commit c73634b; its sole matching code line describes an observational
pneumococcal/J13 edge source (PMID37291502) with2017–2018follow-up,
not NHANES participant-data retrieval, cohort-building or prediction results. Existing beta-carotene literature mentions in
data/diagnostic_panel remain as already disclosed above.

For reproducibility, the exact sorted98refname/SHA pairs are stored below.

```text
refs/remotes/origin/agent/publication-companion-redesign d4b4747fe6eb706d21b12971e72dabc3b9804527
refs/remotes/origin/agent/separate-furg-aptamer-mortality-registries 061506dc995316c6638e8dd4d767c247e58d6d3a
refs/remotes/origin/automation/living-evidence/34119100636-1 3aa4bbfb5f62fd2430686a628ed1f21f511ebd15
refs/remotes/origin/automation/living-evidence/34841978232-1 b2733805abddee01b70f798d49b5852c762b4eca
refs/remotes/origin/automation/living-evidence/36427324812-1 0bfead2d8c8e7c9573d18c8294297a6052abe911
refs/remotes/origin/automation/living-evidence/37320964058-1 fd6f8433f3bf65aaec88568f00fbd5ddc5412693
refs/remotes/origin/automation/living-evidence/38265431170-1 9cc7af6e128aecf30824e017514fd783d3225363
refs/remotes/origin/automation/living-evidence/39252776450-1 68a16e37fad948f50f64f5161e0c69fa96f4613c
refs/remotes/origin/automation/living-evidence/40322874916-1 5369fdb6bdc368c54dbd1935c862da7e321c5373
refs/remotes/origin/claude/absolute-risk-display 87021075c380022bddb2ed4b13c81d36248e5e6e
refs/remotes/origin/claude/aging-biomarkers-diagnostic-system-2i3128 465c235ef9f5807b61f035e86dee396468514e94
refs/remotes/origin/claude/aging-biomarkers-review-publish-g5efzg 99ab87babb6af255591d68a048ca27ea5972bf0d
refs/remotes/origin/claude/aptamer-conference-slide-ivdrv4 f25b64fdf697041d1ab38b3f4dd6076186708a5f
refs/remotes/origin/claude/biomarker-age-assessment-review-s41ya1 f819eb6be8eea7ac50ed1b55352a8c7e26378149
refs/remotes/origin/claude/busy-gates-gm1grs 83a45a8b1b6cb0632e24bdf22d3ed5150deb29f2
refs/remotes/origin/claude/cholesterol-mortality-paradox-mmdv0w 72b97f837c60c6f11785baf1c55b8a89984b7904
refs/remotes/origin/claude/cross-country-mortality-data 0ab03bc4e21b9830c527dd719e60eceb9010ddae
refs/remotes/origin/claude/expert-panel-iteration-10 7421eb8120ce3ab455b1aed3f25df2f9f6e383b6
refs/remotes/origin/claude/expert-panel-iteration-11 27ebaf18a0b22c0995e2d46ac234e41e5caf5fe0
refs/remotes/origin/claude/expert-panel-iteration-12 de77aec47abcb7010261c83353d4b9b1934aa9cb
refs/remotes/origin/claude/expert-panel-iteration-13 22a53f978f1077753d7fab6dcd76d8b19198c711
refs/remotes/origin/claude/expert-panel-iteration-14 36cf93daa803ce044b7bf8c6ad39527e019e41ab
refs/remotes/origin/claude/expert-panel-iteration-15 02c0ed84b255437b192b364096729f9cb100f586
refs/remotes/origin/claude/expert-panel-iteration-16 ea23921f27ee6b952ac47c0147ef3bdc376e9101
refs/remotes/origin/claude/expert-panel-iteration-17 48af02cf3ecec238fe66667b8b46ca7cb52b8967
refs/remotes/origin/claude/expert-panel-iteration-18 f18b3e63b57ff4436eb9ef1c9cf220be4bf826cc
refs/remotes/origin/claude/expert-panel-iteration-19 fa1b89cc996e1fc6d586dafe7c13d093eedc7cac
refs/remotes/origin/claude/expert-panel-iteration-2 633844e2213604d933efe6e5e2f71dbe537c1d7b
refs/remotes/origin/claude/expert-panel-iteration-20 0a00bd461cbb25b6cb70150d52af22e64f40b194
refs/remotes/origin/claude/expert-panel-iteration-21 7a489de788cec433432d4297ae8467d850d5d5dc
refs/remotes/origin/claude/expert-panel-iteration-23 5d83d9fc046a94bed70a01ee3509652456a1fb5b
refs/remotes/origin/claude/expert-panel-iteration-23b 0fde33047037b7a20fb8e962dfccda31954af035
refs/remotes/origin/claude/expert-panel-iteration-24 4c0cb3e182129e09c79beba00da11b2666ad0a6b
refs/remotes/origin/claude/expert-panel-iteration-25 22a1cf50545395d842d28d083c18d2ce08b8856c
refs/remotes/origin/claude/expert-panel-iteration-26 c573d84757cdb97ded4e10bf94fd60cc507726b3
refs/remotes/origin/claude/expert-panel-iteration-26b 35a6cb04b5bf4af777d7aa8163cfd94b2b230bda
refs/remotes/origin/claude/expert-panel-iteration-27 415c3694a05ae74cca68c0e60586b1124ffa39ad
refs/remotes/origin/claude/expert-panel-iteration-28 84df9b105dab5d705876cbe7703fcc216eb731dd
refs/remotes/origin/claude/expert-panel-iteration-29 eb789ee9adcd783de0897ba6c1904602cc7c1f5f
refs/remotes/origin/claude/expert-panel-iteration-3 b4ff8d6345c15e1addd174addb1656882d29441e
refs/remotes/origin/claude/expert-panel-iteration-30 2a85f7661e6c5a9ee8163c1473cd8cfdacd654ed
refs/remotes/origin/claude/expert-panel-iteration-31 71e28fd213938d73ac5c034615eb85ebf99aaee8
refs/remotes/origin/claude/expert-panel-iteration-32 bb8a57adbbdf1e9517beb8036cb46151bdb834de
refs/remotes/origin/claude/expert-panel-iteration-33 81fe1cfafdf4c3345072aa2326579abe22dfae8c
refs/remotes/origin/claude/expert-panel-iteration-34 a6d389ad4043294cd9b09abd3394d3b587aed64a
refs/remotes/origin/claude/expert-panel-iteration-35 f49018f693bbdab46dd4c3586eb42d0b9bbe5efe
refs/remotes/origin/claude/expert-panel-iteration-36 d6404828e92a3d7bb3fda3bee1e9154a59b34819
refs/remotes/origin/claude/expert-panel-iteration-37 5ad2ad635f820e249b396f16649dedf20909deb5
refs/remotes/origin/claude/expert-panel-iteration-38 302029c00c71abfa021e9297aa83cf38fea51361
refs/remotes/origin/claude/expert-panel-iteration-39 19215b31c5f1b8dc540d65ee1b194ec62af7236f
refs/remotes/origin/claude/expert-panel-iteration-4 c417824e0e103c0ddd5aaefdf4622933b747cd5a
refs/remotes/origin/claude/expert-panel-iteration-40 b99c08f1b274c79b77c6a7625cf9a9d623014966
refs/remotes/origin/claude/expert-panel-iteration-41 e72da1929a50beb04598059390b7704bbd484a63
refs/remotes/origin/claude/expert-panel-iteration-42 bfa88478f60cc4304f3522f23d263a30dc967787
refs/remotes/origin/claude/expert-panel-iteration-43 a8e4dddc88a1f559ee2ff50cd0d9c3bdb7cd9921
refs/remotes/origin/claude/expert-panel-iteration-44 7e449aff9496442972228686278ac9240b9e3c9c
refs/remotes/origin/claude/expert-panel-iteration-45 4266b537a76645af7ae93a9f71dc0539500b0c4d
refs/remotes/origin/claude/expert-panel-iteration-46 465fbc482036140f2280a772767c838a85f768d2
refs/remotes/origin/claude/expert-panel-iteration-47 46cf3b6200e9a553843d6d24d31bfa491153969a
refs/remotes/origin/claude/expert-panel-iteration-5 edde1703fcfb1473ee46b63b4fc0c3c52c27cc64
refs/remotes/origin/claude/expert-panel-iteration-6 3c9d3eade93a550cf5b74128300b5ca0b9d89556
refs/remotes/origin/claude/expert-panel-iteration-7 1b249e6e02d5a00d61420668e49f73b0d67fe5ef
refs/remotes/origin/claude/expert-panel-iteration-8 6f6687c8f6ddcace9b95109b9616ac52ca980b6c
refs/remotes/origin/claude/expert-panel-iteration-9 92f32e076643463245c1fad0119eef011fdedd12
refs/remotes/origin/claude/fix-crp-duplicate-marker 51797cc1d3f47c5c05d0a7c4f4740ef5e905053d
refs/remotes/origin/claude/fix-pdf-parsing-crp-wbc cfbc25b905bdd174cab98d0561945df0f037e4c3
refs/remotes/origin/claude/general-risk-levers 9a6aabfd383a142e21fbf7ed300d74416d5bf91d
refs/remotes/origin/claude/gompertz-country-chart 4b047201bc00166c6583545ab6df19bbd2c64a92
refs/remotes/origin/claude/kodeks-login-graph-nrji4b 362021ac49d36b0369e8f9378e7bde2947cca92b
refs/remotes/origin/claude/land-general-levers 59309be9324529552613e750379d89041dcaa485
refs/remotes/origin/claude/mortality-risk-analysis-8o4j37 9a7f236ca7bfaa092e283e329b2796159e2d6b25
refs/remotes/origin/claude/nature-science-article-7thk1y 352088d71fe7449cf5f1445f31f768b47c1c33c6
refs/remotes/origin/claude/pdf-history-multi-upload 2fe2c31bec92b6c3fe1dde0ae3185735cbc7347f
refs/remotes/origin/codex/bidirectional-diagnostic-panel-2026-09-08 3e4ea3a5cb150526cfbab2ab4cfd3fe60a188736
refs/remotes/origin/codex/conference-launch-audit-2026-09-09 dcbf0a8775b2f9b92090420c1be45aeed6615d2f
refs/remotes/origin/codex/conference-run-brief-2026-09-09 f37b3fbda2373c714f1907575ff494d1bbcc4132
refs/remotes/origin/codex/evidence-model-integration-2026-09-15 1a6bd86d125db62dcf5ae632161166c49bf83151
refs/remotes/origin/codex/integrated-theory-audit-2026-09-06 cca645076431fd4845d3e0f323f643f46ed5a109
refs/remotes/origin/codex/iygf-abstract-revision-2026-09-07 4aeb8d96578d4b5f2c6c459cc38f151398c00a41
refs/remotes/origin/codex/laboratory-panel-delivery-2026-09-08 10fde02e67268522f7d8d04e45b49846118a4b57
refs/remotes/origin/codex/living-evidence-theories-2026-09-06 d36d0d1bf2c04f0a205eb4271fb60e39a1a858b2
refs/remotes/origin/codex/panel-depth-audit-2026-09-08 f7502e7be1f53b7a50696af61b66ba079dace4e0
refs/remotes/origin/codex/panel-evidence-completion-2026-09-09 4b491fc8b6fa62d78863be1be05c06384137eb69
refs/remotes/origin/codex/panel-run-reliability-2026-09-08 118147d2c4a57031749ac3174404b13f056fe71d
refs/remotes/origin/fix/pdf-unit-autodetection 010911b524e1d7bd2afa3e2bd25904c50f578628
refs/remotes/origin/fix/result-overlay-layout 194e31dc2f78dc1a7038e06b083e995166c052a4
refs/remotes/origin/fix/unit-audit-and-bp 5a628938d16fbabb3dfa573174491e4896d5a53b
refs/remotes/origin/main 51b3670524951af080634c42ed172f935f17d0cb
refs/remotes/origin/panel-depth 0a65bbb17514dfaaabb1881b7f965028a3a2024f
refs/remotes/origin/panel-laboratory d9d5fa0abfe11446acde2f21370566fd4d8cffa6
refs/remotes/origin/panel-next-run 2a612eadd15289006f76a3d413c3f84ccf132197
refs/remotes/origin/research/panel25-status 801e85972a648951c9922d7715a7728d04d0c11f
refs/remotes/origin/run/conference-panel-2026-09-09 0ea9ad2807086d6d515af5b63ea00ca5ba58cf80
refs/remotes/origin/ux/anatomical-organ-map 35ea005284a79752a7fd4433bb7b6ba3798cb2ad
refs/remotes/origin/ux/disease-card-filters 7b1c47ac0e92687c5faaf670d81dc0f4e9bff1ac
refs/remotes/origin/ux/pdf-review 971076090d84beb183c1672435bf4823db01cc4c
refs/remotes/origin/ux/risk-language a0e7e27edf4079803448f8c3fdbed8688083a238
refs/remotes/origin/ux/trust-boundary-v1 7719bf4764762beca0ed9b73b3cd4a07dee00a0c
```
