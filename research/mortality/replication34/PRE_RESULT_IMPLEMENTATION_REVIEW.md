# Pre-production corrections and numerical validation contract

No production scenarios/MonteCarlo cases, new participant outcomes or predictive
fits were executed before this review. Protocol4df636e and its12scenarios/4laws,
seeds, estimands and neutral33B result are unchanged.

Independent mathematical review found two real boundary problems before
production. Fixed absolute1e-14 coverage slack could force100% coverage for a
legitimate tiny-variance law whose exact normal-z coverage is87.5%. Positive-
variance coverage now has no such slack; only a truly zero-variance law has an
explicitly scoped, scale-relative constancy check. Exact tiny-law and degenerate
fixtures test both branches. A second absolute moment-feasibility tolerance
could clip positive variance tozero at mean±1. The stable bound
(1-|mu|)(1+|mu|), relative8machine-epsilon roundoff only, and positive-variance
underflow rejection fix this without changing scenario inputs.

The already-computed fourth centered moment of simulated sample means is
exported as one aggregate scalar so the independent auditor can recompute the
approximate MonteCarlo variance-SE. No simulation arrays are exported.

The runner's numerical agreement checks use a declared6MonteCarlo-SE band for
meanD,varianceDbar andmean event count. These check stipulated-law arithmetic,
not statistical power or real-data model validity. Coverage is reported with
MonteCarlo uncertainty and is never forced to95%; the rare-event and degenerate
cases can legitimately have different coverage. Input cases/seeds will not be
changed in response to these results.

Independent source review also clarified three metadata descriptions before
production. UKB40001 establishes an underlying ICD10 field, not an already
frozen source-specific crosswalk/missing-cause policy: both mapping checks are
not_verified. NHANESIII repeat's absent population weights alone do not disallow
a justified received-domain refit; sampling contract is not_verified while
programme reuse remains verified_false for independence. UCI's public schema
was previously reviewed, not imported/fitted: remove the erroneous generic
'used by programme' wording. None changes the zero-qualified-source conclusion.
