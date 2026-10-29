# Final implementation freeze before empirical counts

The protocol remains82eacda12eced3c0105019fee09236b647821f45, unchanged.
The initial implementation3b0ca8d was reviewed and all28 synthetic tests passed.
Before opening the empirical PR, a small final patch additionally checks that
the groups object is a dictionary and its ten canonical cause IDs match the
retained definitions. The same existing synthetic test now also rejects an
all-CVD replacement for heart and a missing groups object. No domain, endpoint,
weight, ranking, support threshold, partition or computation changed.

This patch is retained as a separate pre-count commit. No empirical gate or
new cause-specific counts, fits or predictive metrics were executed or viewed
before it. The earlier review's gate-file SHA describes its reviewed snapshot;
the empirical output will export the exact final99 scientific file hashes.
The28 tests were repeated because implementation changed. Joint integration
and Fine–Gray remain deferred; all previous95 scientific modules are preserved.
