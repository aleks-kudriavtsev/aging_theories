# Source interpretation addendum before production

The frozen access/assay reviews are retained as evidence of the review process,
not quietly rewritten. The executable18-release candidate ledger incorporates
the following stricter clarifications before any production qualification or
scenario run:

- FHS: the measured-protein lead is **Offspring exam7**. The Omni1 visit number
  is unverified and is not asserted to be exam7. The earlier review's combined
  Offspring/Omni label is not an approved same-visit package. Current concentration
  units/QC,exact routine11 and underlying-cause mappings remain unqualified.
- UKB: field40001 is underlying ICD10, but current package/scope/code-cleaning/
  crosswalk are not frozen. Both target-mapping states are not_verified.
- NHANESIII repeat: a convenience sample without population weights could have
  a declared conditional estimand; absent weights alone are not a categorical
  failure of algorithm replication. Its sampling contract is not_verified;
  participant/programme reuse remains the independent-source hard blocker.
- UCI: public download/schema documentation does not mean records were imported
  or fitted by this programme. Its protein/mortality endpoint failures remain.

These are metadata/protocol clarifications, not new outcome findings or changed
criteria to admit a favourable source. No access request/agreement/account or
new participant download occurred. The exact updated metadata bytes are pinned
in SOURCE_LOCK.json and will be audited against the scientific commit.
