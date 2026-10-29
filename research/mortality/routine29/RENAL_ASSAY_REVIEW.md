# Primary renal/routine metadata review29 — 2026-10-29

Metadata-only review. No person-level files, source outcomes, model comparisons,
accounts or agreements were accessed. No repository files changed.

## Sources actually read

- Official SAS layout: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/lab.sas
- Official codebook (lowercase URL required): https://wwwn.cdc.gov/nchs/data/nhanes3/1a/lab-acc.pdf
  Catalog76300, version1.1, revisedSeptember2006,135PDFpages.
- Calibration study author-institution record/abstract:
  https://pure.johnshopkins.edu/en/publications/calibration-of-serum-creatinine-in-the-national-health-and-nutrit-4/
  Selvin etal2007, AmJKidneyDis50:918–926,
  doi:10.1053/j.ajkd.2007.08.020.
  Publisher full text returned403; no bypass attempted. The exact equation is
  independently printed in the official CDC codebook, so no secondary article
  or inaccessible table is needed to specify it.

## Creatinine: exact direction, units and limits

Official CEP note, zero-based PDFpages110–111 (textlines4924–4952), explicitly
recommends correction of original NHANESIII serum creatinine.190stored specimens
were compared. Original measurements used kinetic Jaffe alkaline-picrate on
RocheHitachi737. Comparator was CCF Rochecoupled enzymatic assay on RochePModule;
calibrators were traceable to IDMS/NISTSRM967 and checked using CAPLN24.

The official Deming equation is:

`standard_creatinine_mg_dL = 0.960 * original_CEP_mg_dL - 0.184`.

Do not reverse it, apply it to urine creatinine, apply it twice, or call the
original CEP already standardized. Apply after missing-code removal. The
reported CEP frequency range is0.3–13.9mg/dL, not an imposed clinical validity
range. Any nonpositive corrected result must remain missing/stop log processing,
not be floored to a small invented positive value.

This supports a variable explicitly called *CDC-equation calibrated serum
creatinine*, with a primary log2 transform if predeclared. It is NOT measured
GFR and need not be converted to eGFR. A new within-cohort core can avoid an
eGFR equation entirely, keeping age/sex already in the clinical core. No eGFR
assumptions or duplicated derived-age/sex renal feature are required.

The author abstract confirms that the published calibration specimens were
approximately200 per survey, age>=60, with limitations for within-survey drift
and freeze/thaw. Those limitations remain; the correction is not individual
sample validation or proof of comparability in all subgroups.

## Urine albumin/creatinine: exact units and explicit censoring

Official codebook zero-basedPDFpage95 and URPnotepage122:

- URP is urine creatinine in mg/dL, INPUT1956–1960.
- UBP is urine albumin in ug/mL, INPUT1965–1970.
- `UACR_mg_g = 100 * UBP_ug_mL / URP_mg_dL`.

This is unit algebra:1ug/mL=1mg/L and1mg/dL creatinine=.01g/L. Do not use the
SI urine creatinine field with this100factor or use urine albumin alone as UACR.
Only finite positive appropriately decoded operands may enter the ratio.

Important primary-source distinction:

- URP7.9 is an explicitly published below-detection/statistically-suspect code,
  not an ordinary measured7.9. The note says analytic method LOD1mg/dL but all
  values<10mg/dL were considered statistically suspect and coded below limit.
  Codebook has83such records, lowest ordinary reported value11.3. Do not infer
  a new measured denominator, a new LOD formula or normal hydration from7.9.
- UBP0.4 is a released below-detection replacement (386records), not a measured
  value proven aboveLOD. Appendix2 gives urine-albuminLOD0.5ug/mL; generic
  codebook notes below-LOD replacements asLOD/sqrt2 with released rounding.
  Retain the released0.4 if a reported-values design chooses to include it;
  do not re-impute it, silently claim it measured, or divide bysqrt2 again.
- Primary scientific choice should explicitly address suspect URP. A conservative
  model-ready domain may requireURP>=10 and positive albumin, retainUBP0.4 as
  a flagged reported fill, with a prespecified limited sensitivity excluding
  albumin nondetect fills if desired. Alternatively retain all released values
  under an explicitly reported-values sensitivity, but label the denominator
  uncertainty. This is a protocol decision before outcomes, not data-driven
  range trimming.
- Missing raw sentinels:URP88888,UBP888888, plusblank. Remove before log/ratio.
  Other numeric values must not be turned into missing because of a later
  empirical distribution or clinical risk cutoff.

## Glycaemic marker

GHP INPUT1861–1864 is *glycated hemoglobin(HbA1c),percent*; GHPMETH1865 identifies
Diamat instrument1/2/3 or alternateaffinity4. MissingGHP8888/blank, method8/blank.

The official GHPnote, zero-basedPDFpage113, states that Diamat ion-exchangeHPLC
and the alternate affinity method were BOTH standardized to the DCCT reference.
Affinity was used for problematic chromatograms, variants or degradation;
it is not arbitrary unstandardized total glycohemoglobin. Thus continuous
releasedGHP% is a defensible single glycaemic covariate on the DCCT scale;
GHPMETH can be retained for method/QC descriptive counts. Do not claim the
historical assay is newly NGSP-certified merely from this note. No IFCC
conversion is required or supplied. HbA1c does not equal a complete diagnosis/
treatment history; adding a duplicate diabetes-history flag is unnecessary for
this compact core question unless separately justified.

## Parsing: FORMAT is not INPUT

SAS INPUT uses ordinary column ranges (no implieddecimalinformat). FORMAT CEP6.1,
GHP6.1,URP7.1,UBP8.1 controls display; it is NOT permission to divide raw numbers
or missing sentinels by10. Raw numeric tokens include their decimal point.

SELECTEDINPUT positions:SEQN1–5;CEP1784–1787;GHP1861–1864;GHPMETH1865;
URP1956–1960;UBP1965–1970. FILENAME LRECL1979 includes CRLF; last nominal input
ends1977. INFILELABMISSOVER permits absent trailing unselected content; selected
reader needs all columns through1970. Official release n29314. Decode numeric
literals including decimals and signed/blank syntax explicitly; compare a second
reader before joins. Adult and LAB demographic/design positions differ, so do
not reuse adult offsets for fields other than matchedSEQN.

## Useful predeclared core

Previous verified age/sex/smoking/SBP plus BMI(if separately source-verified),
continuousGHP%, log2CDC-calibratedCEP and log2UACR is compact and interpretable.
No Cys-derived eGFR belongs in the comparator when testing addedCys; it would
mix comparator and exposure. No physician CHD/angina proxy is required. All
comparative models must be re-estimated on the SAME complete-case domain,
using the protein-subsample WTCYPEX6 rather than generalLAB/MEC weights.
