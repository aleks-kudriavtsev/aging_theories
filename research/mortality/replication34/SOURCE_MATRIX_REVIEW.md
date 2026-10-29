# Replication34: preproduction source-matrix review

Reviewed 2026-10-29. Read-only repository/source-evidence review; no network GET,
source records, new outcome inspection, gate execution, SCENARIOS execution,
downloads, accounts, agreements, fits or repository edits. This review is
separate from the programmer's gate-module and precision-code review.

## Scope and snapshot

Read all 18 rows of `replication34/CANDIDATES.json`, all 32 references in
`SOURCE_MANIFEST.json`, `GATE_CONTRACT.json`, frozen ACCESS_REVIEW.md and
ASSAY_REVIEW.md, STUDY_PLAN.md/METHOD_DESIGN.md and the retained freeze30 reviews.
Read gate.py only to understand what a positive evidence state means; did not
import it or invoke its qualification functions.

Initial snapshot SHA256:

- CANDIDATES.json: `aa9c8c7b62797e557875fa2addbb6a80ff5212d06f9df957c39ac77077886b34`
- SOURCE_MANIFEST.json: `31dbf4489edd0a698297d573e16551b8d58948baba016000590be49813ec6370`
- GATE_CONTRACT.json: `affa8b79e61948a5d8da4cd981026aa9f731f1c603f284381f035370d445b541`

Static checks: 18 distinct candidate IDs; 32 distinct source IDs; no unresolved
candidate source references; all rows have exactly the 15 primary plus 2
transport-only check keys. All six review hashes and the prior catalog hash
match local bytes. Original-source byte digests are explicitly null and the
manifest explicitly disclaims having saved or verified original source bytes.
Review-file fingerprints do not imply archived remote-source reproducibility.

## Findings to correct before source production

1. **UKB mapping evidence is ahead of the documented qualification.** The two
   `underlying_*_mapping` states were verified_true (initial lines2085–2097),
   while their reason says the exact study-specific crosswalk must still be
   frozen. The official 40001 metadata establish underlying ICD10 availability,
   distinct from contributing40002. That is a strong endpoint lead, but does
   not establish the exact release-specific mapping, code cleaning, missing
   cause handling and scope. Use not_verified for these two checks until that
   work is fixed. Preserve the positive fact that raw underlying ICD10 exists.
   Parent accepted this refinement before production.

2. **Repeat III design failure was overgeneralized from missing population
   weights.** The initial sampling_design_contract was verified_false because
   it is a convenience repeat sample without population weights (lines244–249).
   In this new algorithm-refit contract a justified conditional received-domain
   estimand/unit-weight design is possible in principle. No population weights
   alone establish neither a valid design nor an impossibility of any valid
   design. Change this check to not_verified, with exact design/estimand still
   unqualified. Its reused-person/source independence check remains
   verified_false, so this does not create an independent cohort. Parent
   accepted this refinement before production.

3. **UCI reason should not imply a previous raw-data model fit.** The public
   access reason was generated as "already documented and used by the
   programme" (initial line2621), but the retained cohort review says records
   were not downloaded. Say publicly documented anonymous release, previously
   assessed as metadata only. The proteins and mortality/time/cause schema are
   explicitly absent in this 24-feature classification release and those
   negative checks are justified. Parent accepted the wording correction.

4. **FHS visit label needs narrowing or explicit uncertainty.** The row says
   "Offspring/OMNI1 exam7", and the human action asks for an exam7 panel. The
   frozen review supports an Offspring exam7 measured-protein lead, but does
   not separately qualify the OMNI1 visit number or its inclusion in the same
   BioLINCC Offspring package. Do not infer a shared exam number across
   generations/cohorts. Safest wording is Offspring exam7; OMNI1 exact visit
   and released panel remain unverified. No new retrieval was attempted to
   fill this gap. Assay reviewer was asked to check its retained evidence.

These findings concern evidence labels and actionable visit definitions. They
do not make an otherwise ready source fail or turn any source ready. Each row
already has an independent unresolved or false primary requirement. Expected
ready count is therefore zero by inspection; production gate was not executed.

## Material boundaries correctly represented

- The new primary mode is algorithm replication with refitting in a distinct
  source, not transport of a single33B frozen coefficient vector. The40OOF fit
  records are not one preselected final scorer. Transport-only artifact/unit
  checks remain separate from primary requirements. Algorithm-level precision
  is explicitly unqualified by a conditional fixed-loss planner.
- BTP is not required. A different manufacturer of interpretable absolute
  concentration Cys/B2M is not alone a reason to reject a new refit. Matrix,
  units, analyte identity, quantification, calibration/QC and release context
  still need qualification. UKB NPX does not silently become B2M mg/L.
- Public001/002/010 is sufficient in principle for cancer/heart endpoints;
  known-other010 is a competing noncase, not unknown cause. Late NHANES is
  blocked separately by public2019 follow-up, unqualified joint pair and
  programme-use history, not by requiring all ten cause codes.
- Study pages, codebooks and papers are correctly distinguished from actual
  anonymous person-level releases. ARIC/CHS/MESA/FHS are controlled and not
  requested/granted; CRIC access remains unknown after the retained metadata403.
  HRS provider-prohibited AI processing is not described as fixed by a login.
  UKB requests require current availability/scope review; no reopen is assumed.
- ARIC visits2/4/5/6/7 are distinct rows. Missing V5 B2M is unqualified rather
  than globally absent; V2 UACR/V4 HbA1c and late-visit96month gates remain
  unresolved/blocked. MESA2012–14 kidney measurements are not assigned to2000
  baseline; events2020 do not supply every late-baseline survivor eight years.
- CHS consent scope does not automatically allow cancer analysis; linked
  legacy age/BMI processing is a caveat, not asserted to describe the current
  approved package. Late urine/HbA1c substitutions are not authorized.
- Distinct recruitment source and programme-use review are explicitly not
  proof of zero person overlap. No identity linkage/reidentification is proposed.
  Previously used NHANES/III repeat observations are not an untouched cohort.
- Named future human actions require eligible PI/institution, actual approved
  package/operator/environment and source-specific mappings/design. No action
  has been initiated, and no account/DUA is implied to authorize AI processing.

## Review disposition and blockers

**Final read-back: all four findings are resolved in current input bytes.**
UKB both mappings and repeat-III design are now not_verified with appropriately
limited reasons; UCI no longer implies previously analysed records; FHS is
explicitly Offspring exam7 with Omni visit unverified. The future exam7 request
therefore has the stated Offspring scope. No outstanding material source-matrix
classification or human-action contradiction was found after these changes.

Final read-back SHA256:

- CANDIDATES.json: `c6996840da99902d851faca323d60d1ea2a18bdf11c01be1816ffff103090706`
- SOURCE_MANIFEST.json: `08658be210a97b52c642f1b5277fb94af9ac4f5e723efe3ffc961668e2f35b50`
- GATE_CONTRACT.json: `affa8b79e61948a5d8da4cd981026aa9f731f1c603f284381f035370d445b541`

The 18/32 counts, all current review fingerprints and prior catalog fingerprint
still match. No production qualification, comparative empirics or simulation
was executed or authorized by this audit. The underlying empirical blocker remains:
no presently authorized independent joint routine+Cys+B2M+known8year
heart/cancer person-level package is qualified. Relevant controlled requests
remain not_requested/not_granted, with source-specific scientific gaps stated.
No new operational denied/pending/uncertain action occurred during this audit.
No labreportv3, presentation, previous model/result or repository file was edited.
