# NHANES 1999–2004 serum creatinine: cycle-specific qualification

Checked **2026-10-29**, official CDC documentation only. No participant records, outcomes, XPT/raw downloads, accounts, agreements or repository edits.

**Implementation decision:** apply the CDC1999–2000 equation only to that cycle's original serum-creatinine measurement; use the released2001–2002 and2003–2004 serum values without another IDMS calibration. **Do not apply the NHANESIII equation to any of these cycles.**

| Cycle and release | Exact released variable in cycle codebook | Unit | Additional standardization for analysis | Primary justification |
|---|---|---|---|---|
|1999–2000 `LAB18` | `LBXSCR` | mg/dL | `Scr_standard = 1.013 * LBXSCR + 0.147` | CDC strongly recommends correction. Comparison of196stored specimens found significant bias of the original RocheHitachi917 kinetic alkaline-picrate/Jaffé result relative to the CCF Roche enzymatic reference. |
|2001–2002 `L40_B` | **`LBDSCR`**, not `LBXSCR` | mg/dL | `Scr_standard = LBDSCR` (identity; no added IDMS recalibration) | CDC comparison of194stored specimens found the difference within assay error and concluded no correction necessary. Original measurement used BeckmanLX20 and RocheHitachi917 Jaffé methods. |
|2003–2004 `L40_C` | `LBXSCR` | mg/dL | `Scr_standard = LBXSCR` (identity) | CDC comparison of190stored specimens found the difference within assay error and concluded no correction necessary. Original method was BeckmanLX20 kinetic alkaline-picrate/Jaffé. |

The reference direction is **original NHANES serum result → CCF Roche coupled enzymatic result**, with calibrators traceable to isotope-dilution mass spectrometry and NISTSRM967, checked with the CAPLN-24 linearity set. It is not the inverse regression. For a hand check, original1999–2000 serum creatinine1.000mg/dL becomes1.160mg/dL. The same released1.000mg/dL in2001–2002 or2003–2004 remains1.000mg/dL under CDC's no-additional-correction recommendation.

The2001–2002 codebook separately describes a cross-over study harmonizing five laboratory tests, including creatinine, between the two participating laboratories; those regression adjustments were **already applied in the release**. That harmonization and the later no-additional-IDMS-correction decision are distinct. Do not repeat the cross-over adjustment.

All three codebooks also provide `LBDSCRSI` inµmol/L. Prefer the cycle's mg/dL variable above to avoid unit ambiguity. The2003–2004 processing note explicitly defines `LBDSCRSI = 88.4 * LBXSCR`. If aµmol/L source is intentionally used, convert to mg/dL before applying the mg/dL calibration intercept; do not add0.147 toµmol/L values. These are **serum-only** equations and must not alter urinary creatinine or UACR denominators.

## Source discrepancy and required import assertion

The official1999–2000 biochemistry **method PDF** lists `LBDSCR` in its public-release analyte table, whereas the current cycle-specific **LAB18 codebook** lists `LBXSCR`. This is a documentation inconsistency, not evidence that both columns exist. Freeze the exact intended XPT release and assert its actual column against the cycle-specific codebook; stop if it does not match. Do not guess an alias or silently substitute a column. For2001–2002, both the L40_B codebook and creatinine method PDF agree on `LBDSCR`; for2003–2004 both agree on `LBXSCR`.

## Official sources and short quotations

1. [CDC LAB18,1999–2000](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/LAB18.htm), Analytic Notes and `LBXSCR` codebook. Exact equation quotation: “Standard Creatinine (Y) = 1.013*NHANES Creatinine (X) + 0.147”. Direction,196specimens and reference traceability are described immediately above it.
2. [CDC L40_B,2001–2002](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L40_B.htm), Analytic Notes and `LBDSCR` codebook. Quotation: “Correction for Serum Creatinine for NHANES 2001-2002 is not necessary”. The displayed diagnostic regression `CCF≈−0.005+1.010*NHANES` is **not** CDC's recommended extra correction.
3. [CDC L40_C,2003–2004](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L40_C.htm), Analytic Notes, processing and `LBXSCR` codebook. Quotation: “no correction is necessary for serum creatinine values in NHANES 2003-2004”. Likewise, do not apply its diagnostic comparison regression `CCF≈0.02+0.98*NHANES` as a new calibration.
4. [CDC1999–2000 Lab18 method PDF](https://wwwn.cdc.gov/nchs/data/nhanes/public/1999/labmethods/lab18_met_biochemistry_profile.pdf), printedpp1,3,6: Hitachi917/Jaffé method; conflicting public-release field label noted above.
5. [CDC2001–2002 Lab40 creatinine method PDF](https://wwwn.cdc.gov/nchs/data/nhanes/public/2001/labmethods/l40_b_met_creatinine.pdf), printedpp1–2: `LBDSCR`,mg/dL; BeckmanLX20 kinetic alkaline-picrate method; two-laboratory release context.
6. [CDC2003–2004 Lab40 creatinine method PDF](https://wwwn.cdc.gov/nchs/data/nhanes/public/2003/labmethods/l40_c_met_creatinine.pdf), printedpp1–3: `LBXSCR`,mg/dL; BeckmanLX20 kinetic alkaline-picrate method.

Remaining limitation: this review verifies metadata and population calibration guidance, not each specimen's analytical validity or calibration uncertainty. It does not imply that numerical frequency ranges are clinical acceptance limits. The actual XPT columns and hashes have not been inspected here.
