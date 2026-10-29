# Prespecified creatinine-based renal representation

CDC NHANESIII lab-acc.pdf CEP note confirms standard serum creatinine(mg/dL)
=0.960*originalCEP−0.184, based on190stored samples compared with an
IDMS-traceable Roche enzymatic method/NIST SRM967. This equation applies ONLY
to serum CEP, never to urinary URP. Nonpositive recalibrated values are missing.
The empirical calibration has uncertainty not represented by model bootstrap.

NIDDK publishes the adult2021 CKD-EPI creatinine equation for IDMS-traceable
SCr,age,sex:
142*min(SCr/k,1)^alpha*max(SCr/k,1)^(-1.200)*0.9938^age*(1.012 if female),
k=.7female/.9male, alpha=−.241female/−.302male. No race multiplier.
Unit: estimated mL/min/1.73m². It is an estimate, not measured GFR, and the
historical creatinine calibration does not itself validate GFR accuracy in
this cohort or in an individual. No cystatin term is allowed in the routine
comparator; that would contaminate the protein-addition contrast.

Use eGFRcr as the sole serum-creatinine representation: min(eGFR,60)/15 and
max(eGFR−60,0)/30. A fixed60hinge permits different low/high relationships
without searching knots. Do not additionally include raw/calibrated creatinine.
Age and sex remain core risk factors; eGFR is a deterministic transformation
of those and one creatinine measurement, not a separate independent analyte.
No CKD diagnosis or individual renal clinical validity is inferred.

Primary sources reviewed before routine29comparisons:
https://wwwn.cdc.gov/nchs/data/nhanes3/1a/lab-acc.pdf (CEPnote)
https://www.niddk.nih.gov/research-funding/research-programs/kidney-clinical-research-epidemiology/laboratory/glomerular-filtration-rate-equations/adults
