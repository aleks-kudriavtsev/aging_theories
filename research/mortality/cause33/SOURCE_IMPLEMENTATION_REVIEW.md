# Cause33 source/wrapper review (read-only)

Protocol82eacda12eced3c0105019fee09236b647821f45. No GET/retrieval, participant-source reads, fits, predictive metrics or repository edits were performed. Read-only code inspection, pure hash verification and five synthetic intake tests only.

## Result

No critical scientific/code blocker found in the inspected source/wrapper/workflow scope. The new count/endpoint/support gate.py is outside this review and is being implemented/reviewed separately.

## Exact checks

- retrieve.py:19–30: metadata resolver correctly selects the five cause33 local metadata basenames through LOCAL_METADATA and resolves root endpoint_contract11.py and qualified legacy paths against mortality/. Pure check_study verified18 metadata hashes and95 prior scientific Python hashes.
- retrieve.py:33–43: calls preserved b2m31 intake only, then binds protocol full SHA, current cause33 lock SHA, exact base-manifest byte SHA and full per-file fingerprints in the wrapper. No new source URL, mirror, PDF retry or ZIP route is introduced.
- retrieve.py:46–58: verifies schema/protocol/lock/base bytes/per-file fingerprint equality before parsing. This helper alone is only wrapper binding; full source authentication remains the preserved loader invoked next.
- qualify.py:23–27: check_study + verify_wrapper precede b2m31.cohort.load; loaded_base must equal the previously bound base object. Preserved b2m31.cohort.py:17–47 verifies unique safe official URL rows, sizes/content SHA, exact declared file/URL inventory and25 original historical data hashes; actual parsers and independent mortality readers run in that preserved load path.
- b2m31.retrieve.py:77–83,107–133: the two previously denied mortality PDFs are records only and never enter static_urls/get_one; failures prevent successful intake return. Cause33 inherits this behavior without overrides.
- qualify.py:31–46: explicitly labels no models/predictive metrics, deferred joint/Fine–Gray, reviewed dictionary without newly fetched PDF bytes; exports the complete base manifest, wrapper and aggregate inherited source-domain audit.
- qualify.py:47–56 serializes both outputs before directory/file creation and emits exact UTF-8 JSON blocks for independent byte reconstruction. No participant source file is copied to the output directory.
- workflow:8–9 read-only GitHub permission;18–20/33–35 checkout persist-credentials false;40–45 uses only cause33 public intake and qualification;46–51 uploads only /tmp/cause33_results rather than raw /tmp/cause33_cdc. Empirical job is PR-only (line29).

## Verification

python -B -m unittest discover -s tests -p test_cause33_intake.py -v:5/5 passed (0.027s). All retrieval calls in these tests are mocked; fixtures are synthetic temporary files. Protocol commit82eacda was inspected as a metadata-only commit with six new protocol/review/lock files, before the new execution module.

The retained dictionary/public endpoint review is intentionally not a new verification of the two previously denied PDF bytes. This limitation is explicit and does not become a fresh source validation claim.

## Reviewed file hashes

- research/mortality/cause33/retrieve.py: fe6da2af9e3295356f11e5b828316ee7770c26b6fbcfffacac7963f3ca5d376d
- research/mortality/cause33/qualify.py: 347908f5705fe943e5e8bbbdee2f06d765051b3f172dd78d317752335af80c2f
- research/mortality/cause33/LOCK.json: 751d175b648f2fedf2b6562a68c779c057e1d2b8035fbd9c8e2eed9b2b18c02c
- .github/workflows/cause33-gate-checks.yml: 1fae8bad5a11ffa5d5fa447bec6fc09f6ae3e512c86ebe1f32329eb15ce92048
