# Metadata-triggered amendment, before new outcome inspection

2026-09-10. Primary plan55279f23 and primary model SHA40069d1a98a29308bee26e5dac895388f4c9008b983bc834ad8ac0a57fbcd094 remain unchanged. The2011-2014 files have been retrieved, but only HTML measurement documentation has been examined at this point, not individual XPORT or mortality records.

CDC CBC_H Analytical Notes report a shift from RDW12.84 to13.55% and WBC7.03 to7.41thousand/uL between2011-2012 and2013-2014, coinciding with hematology instrument changes. This cannot justify a patient-level correction using population means: those are not paired specimen calibration data. Source: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/CBC_H.htm .

Add ONE secondary sensitivity panel omitting rdw and log_wbc from the primary22-feature no-CRP panel. It retains13clinical and7engineered laboratory features. Development and comparison use exactly the same15280 complete no-CRP profiles and the same complete-profile evaluation domains. This avoids changing population selection together with panel composition. No claim to validity on missing-CBC profiles will follow from this sensitivity.

Fit only on1999-2010 with the previously locked3calendar-block folds, penalty grid1/10/100, equal mean heldout5year six-stateBrier and tie rule. Freeze the sensitivity model checksum before individual evaluation inspection. Compare its fixed predictions with both primary panels in2011five-year and2013four-year domains. Use1000 paired PSU draws. Do not replace the primary result with the best of three panels, and do not recalibrate based on evaluation results.

Implementation clarification: independent cycle bootstrap streams use seed20260912+(survey_start-2011), yielding20260912 and20260914. Missing-direction AUC repetitions are counted explicitly. This changes no fitting or hypothesis decision.

This is a secondary analysis added after retrieval and metadata review but before evaluation outcomes, not part of the original pre-retrieval lock. The primary model/code remain unchanged. Exact data/code/model checksums and both favourable and unfavourable results are retained.
