# Cause33 pre-count code review —2026-10-29

Protocol82eacda12eced3c0105019fee09236b647821f45 was fixed before new cause
counts. Reviewed current gate.py, retrieve.py, qualify.py, support/definition
metadata, workflow and synthetic tests. No source records, actual source GET,
comparative metrics or model fits were read or run; no repository edits.

Disposition: **no remaining code/scientific blocker found before empirical gate.**

One discovered issue is FIXED: gate.later_comparisons_allowed originally depended
only on absence of unknown early causes and could be true with no supported
priority cause. Current gate.py279–285 computes proceed_codes first, separates
source_endpoint_qualified, sets later_comparisons_allowed=bool(proceed_codes),
and exposes no_supported_priority_cause when relevant. Boundary test verifies it.

Verified fidelity:
- exact full31 common flag00 domain, native W4/W2 provenance and pooled W6 formula;
- valid original1999/2001/2003 source years, underlying cause1..10, known vital
  status/time,96month inclusive horizon and no short-observed survivors;
- unknown common early deaths remain label-unknown and block every comparative
  path rather than becoming negatives; later unknown deaths are known8y noncases;
- all9 defined causes ranked by weighted event burden/raw events/code, top5
  priority fixed BEFORE support, no backfill;010 accounting never clinical rank;
- ten groups+unknown raw/weighted accounting, identical target denominator;
- fixed31 outer/inner hash maps constructed from full eligible masked units,
  same eligible/common inventory and actual3-PSU stratum;40 support cells and83
  threshold checks per cause, without data-dependent rehash;
- pure gate imports only hashlib/json/numpy/pandas and has no fit/predictive-metric
  imports or calls. Runner exports only aggregate JSON and full source wrapper.

Independent source/wrapper review: /tmp/cause33_source_review.md. Resolver puts
five local metadata names under cause33 and retained endpoints/old files under
mortality. It checks18 metadata pins and95 previous scientific hashes. Source
wrapper binds protocol+lock+exact base-manifest bytes+all per-file fingerprints;
the preserved31 loader then validates historical data SHA/official URLs/readers.
This describes checked code paths, not an independent raw-source audit here.
Workflow publishes only /tmp/cause33_results, not source data, and inherits the
two prior-denied mortality-PDF statuses without another GET.

Verification run: python -m unittest discover -s tests -p 'test_cause33*.py' -q
Result:28 synthetic tests,0.690seconds,OK. Includes source wrappers, runner
atomic serialization/changed-manifest checks, unknown causes, competing noncases,
rank/support boundaries, folds and no-fits/no-participant-export behavior.
Empirical data qualification is still pending; no cause count/ranking claim is
made by this review. Later comparative fitting still needs its separate committed
pre-comparison plan; user research authorization is already present. Joint
integration and Fine–Gray remain deferred.
