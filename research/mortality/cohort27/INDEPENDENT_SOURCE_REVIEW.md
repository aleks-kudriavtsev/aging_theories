# Cohort27 source/endpoint/assay implementation review

Read-only review of gate.py, retrieve.py, importer.py, qualify.py, CONTRACT.json, LOCK.json and ANALYSIS_PLAN.md on2026-10-29. No participant records downloaded, no source/model/results changed.

## Verdict

No substantial scientific or technical defect found in the declared numeric-availability/endpoint-count qualification path. Full frozen-model evaluation remains explicitly blocked; no AUC/Brier/risk probabilities are calculated.

## Confirmed against official metadata

- Primary first-exam routes27a/SSCYSTAT and36a/SSNH3BTP are distinct from37a/SSNH3CYS second-exam convenience data. Dynamic exact-directory+basename resolution cannot quietly substitute the repeat file.
- WTCYPEX6 is the published first-exam serum-cystatin MEC subsample weight, reused for B2M/BTP. The importer verifies overlapping observed weights and takes only positive finite weights for its protein-domain qualification. This is not the ordinary MEC weight and not WTCYPHX6(MEC+home).
- Adult selected fields and inclusive positions match official adult.sas:SEQN1–5,HSSEX15,HSAGEIR18–19,HSAGEU20,SDPPSU6 43,SDPSTRA6 44–45. Record content3346characters excludesCRLF. HSAGEU=2 denotesyears; interview age40–75 is labelled correctly.
- Official adult PDF p72 confirms combined pseudo-strata01–49 and PSU pair1/2. The current PSU check is therefore appropriate.
- The original cystatin field and first-exam B2M/BTP field labels are mg/L. Conflicting prose assay-range units are not used for conversion or QC clipping. Positive numeric presence is not promoted to LOQ/PASS; original cystatin calibration is labelled.
- Official public-use LMF description p3 explicitly states that ALL surveys use the SAME file layout. Thus the inherited continuous-NHANES position specification is also applicable to NHANESIII, rather than being an unverified assumption. Description p1 and dictionary p2 confirm PERMTH_EXM originates at MEC/exam for bothIII andcontinuous.
- Deaths through96months are counted; alive records observed<96months are reported as early censored. Known status at8years includes later deaths correctly. No artificial survival label is assigned to early censors.
- Source and metadata gates run before table parsing, source bytes/size/SHA and exact release-index hrefs are checked, selected adult fields are independently parsed bymanual andpandasFWF readers. One-to-one joins enforce identifier uniqueness.

Official metadata URLs:
https://wwwn.cdc.gov/nchs/data/nhanes3/1a/adult.sas
https://wwwn.cdc.gov/nchs/data/nhanes3/1a/ADULT-acc.pdf
https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.pdf
https://wwwn.cdc.gov/nchs/data/nhanes3/36a/SSNH3BTP.htm
https://wwwn.cdc.gov/nchs/data/nhanes3/37a/SSNH3CYS.htm
https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf
https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf

## Nonblocking reporting limits

1. gate.audit returns participant_records_processed=false as a metadata-stage snapshot. When embedded in a successful real qualification result, explain that the value describes the PRE-READ gate; it is not a claim that the entire qualification run processed no participant records.
2. Shared weight checking ignores missing pairs by design and allows the larger BTP weight domain even where cystatin measurement is absent. The reported eligible domain should therefore be named positive published first-exam protein-subsample weight domain; the joint-three numeric subset is the appropriate measured intersection.
3. No standardization/LOQ adjudication or complete frozen clinical-predictor mapping is established. Qualification may count data and events; it cannot unlock predictions merely because n/events are adequate.
4. Source/reader consistency checks do not remove possible linkage-selection bias, stored-serum degradation, original-assay drift or synthetic perturbations of public followup time. Counts are observed published-data counts, not complete long-term-risk validation.

## Reader-width correction after first-run stop, before counts

The first real run reportedly stopped at the strict adult-record width check. My earlier wording that3346was the content width of every physical record was too strong:3346is the nominal final INPUT endpoint, while the official SAS program declares LRECL3348includingCRLF and uses MISSOVER. Those statements do not establish that all distributed physical lines retain trailing blanks through column3346.

The selected six-field contract only requires support through column45=max(end for selected official INPUT positions). A reader can defensibly require45<=physical_line_length<=3346, preserving every character in that selected range, without requiring unused trailing columns. This is an input-format correction, not a change to age, biomarker or outcome definitions. Before implementing it, the runner should report aggregate line-length frequencies and terminal-control/blank classifications, with no raw lines/IDs. Do not silently skip arbitrary short records, truncate nonblank data or broaden the maximum to an observed number without metadata support. An isolated terminalDOS EOF0x1A, if demonstrated by aggregate classifier, is a separate text-terminator case rather than a participant record.

No person records were retrieved or inspected by this review. The actual physical-width distribution and reason for the first-run mismatch remain unknown to this reviewer at the time of this correction.
