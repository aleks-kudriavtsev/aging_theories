# Rep34: scenario precision/event planner, method design only

2026-10-29. No repository edit, empirical model fit, split search, new outcome
read or33B effect-based planning has been performed. This is a design for a
future frozen method/protocol, not an external-validation result or a claim
that a cohort/sample size is clinically adequate. No source coefficient artifact
is being fitted or required merely to run this scenario calculator.

## Recommendation and conditioning

Make the first implementation a **conditional paired-loss precision planner**.
Its analytical N answers only: under explicitly supplied variance/design
assumptions for fixed prediction rules, what evaluable N gives an approximate
normal CI half-width h for a paired difference in Brier loss?

It does NOT answer sample-size adequacy for nested refitting/algorithm
replication. Shared training, penalty selection, preprocessing and overlapping
training folds introduce dependence/instability not represented by an iid
individual-loss SD divided by sqrt(N). For algorithm_replication mode return
`conditional_proxy_only` / `algorithm_precision_unqualified`, not an adequate-N
or power success flag. A later full-pipeline generative simulation would need
predictors, outcomes, clustering and the complete training algorithm; direct
loss simulation cannot validate training uncertainty. If a user explicitly
requests a training variance multiplier>=1, allow it only as a named heuristic
sensitivity, never infer it from33B or present it as a proven correction.

A frozen external METHOD protocol can be prepared now. Distinguish:
1. **Algorithm replication:** prespecified nested training/validation in a new
   cohort tests the algorithm in that setting; it is not transport of a single
   old coefficient vector.33B produced outer-fold fits, not one universally
   frozen source scorer. This planner alone does not establish replication N.
2. **Locked-source coefficient transport:** before accessing external outcomes,
   a separate source-training plan must select/tune on source data only, fit
   and freeze the intended eight cause/arm scorers, transformations, assay/unit
   mappings and hashes. Then evaluate those unchanged paired scorers externally.
   Fixed-scorer precision is conditional on those frozen parameters; source
   training/selection uncertainty is excluded unless separately modelled.
3. Updating/calibrating with external labels is model updating; it needs its own
   training/evaluation separation and cannot be sold as unchanged-scorer
   validation. Diagnostic calibration must not feed back into scored risks.

## Paired estimand and analytical calculation

For each prespecified cause and contrast, take NEW minus REFERENCE:
D_i=(Y_i-p_new,i)^2-(Y_i-p_ref,i)^2, with Y binary and predictions in[0,1].
Negative mu=E(D) favours the new rule. D lies in[-1,1], and
D=(p_new-p_ref)*(p_new+p_ref-2Y).

Let v=Var(D), n=number of jointly evaluable paired participants, and DEFF_D>0
be a user-specified TOTAL variance ratio relative to iid/unweighted sampling:
Var(Delta_hat) approximately DEFF_D*v/n,
SE approximately sqrt(DEFF_D*v/n),
h approximately z_(1-alpha/2)*sqrt(DEFF_D*v/n),
n_required=ceil(z_(1-alpha/2)^2*DEFF_D*v/h_target^2).

Here h is HALF-width; full width is2h. With fixed v the required precision N
does not depend on mu. Do not turn an observed delta or observed significance
into a post-hoc power number. At95% use the exact numerical normal quantile,
not a rounded constant if testing integer ceiling boundaries. This is a normal
approximation, not guaranteed finite-sample coverage, especially for skewed
losses, few independent clusters or unstable weighting. A large participant N
in a few clusters does not cure few-cluster uncertainty; report the planned
cluster/design support separately.

If the input is paired SD s_D, use v=s_D^2. If the input is an assumed/reference
SE, REQUIRE its n_ref and DEFF_ref: v=SE_ref^2*n_ref/DEFF_ref, with explicit
conditioning and provenance. Never confuse an SE of the mean with person-level
loss SD. Separate SDs/SEs for the two Brier scores are insufficient without
within-person loss covariance: v=v_new+v_ref-2Cov(L_new,L_ref).

Validate moment feasibility: mu in[-1,1],0<=v<=1-mu^2 (allow only declared
floating tolerance). Incidence pi ALONE does not determine v. Even at rare pi,
prediction-pair heterogeneity can make Var(D) large; do not impose
v<=4*pi*(1-pi). That expression applies only to a particular conditional-Y
component with fixed predictions, not the whole participant distribution.
For zero v, emit a degenerate illustrative scenario flag; numerical n>=1 is not
proof that either model, calibration or AUC is valid.

Unequal survey weights need an unambiguous DEFF basis. If fixed weights and iid
D are independent, Var(weighted_meanD)=v*sum(w^2)/(sumw)^2=v/n_Kish. Do not apply
Kish inflation AND a total DEFF that already includes unequal weights. Neither
that identity nor a generic DEFF is automatically valid for informative weights,
clustered outcomes or unknown loss/weight dependence.

## Scenario inputs and uncertainty ranges

Use a frozen, versioned scenario file, separately identified by cause and
primary/secondary contrast, with:
- horizon8years and population/endpoint definitions, validation_mode;
- assumed evaluable-domain incidence pi, assumed mu_D (optional for precision),
  paired SD/variance OR fully anchored SE_ref, and their provenance;
- DEFF_D and its definition/range; event-count DEFF_Y separately if used;
- stakeholder/user-selected h_target and alpha; label all supplied numerical
  examples `illustrative_technical_not_clinical_acceptability`;
- JOINT probability of available qualified core plus BOTH proteins, q_joint;
  q_known=P(required endpoint known|jointly eligible), source/sampling/assay
  qualification assumptions, right-censoring proportion c;
- incidence/SD/DEFF/availability LOW–CENTRAL–HIGH scenarios, treated as specified
  sensitivity assumptions, not estimated CIs, priors or probabilities of truth.

Do not automatically import33B deltas, bootstrap SEs or ranked-marker results.
No silent defaults implying that the current historical effect transports.
An independent future pilot/literature estimate can be an explicitly supplied
assumption with provenance/uncertainty and a separately qualified population.
For cancer/heart incidence scenarios ensure mutually exclusive cause incidence
assumptions sum<=1 when presented as one joint population. Do not invent a
cross-cause covariance or claim a between-cause effect difference.

A toy example may use mu_D=0,SD_D=.05,DEFF_D=1.5,h=.002,pi=.03,q_joint=.75,
q_known=1,c=0: it is only an arithmetic fixture, not a recommended clinical
precision target, expected biomarker benefit or derived external event rate.
Prefer an empty user scenario template plus clearly marked test fixtures to a
large grid of apparently authoritative sample-size recommendations.

If every one of four cause/contrast pairs must meet its own declared marginal
half-width on the same domain, take the maximum conditional n across scenarios
for that particular requirement. This does not make their CIs simultaneous.
Familywise targets/multiplicity would need a separately explicit protocol.
Brier precision does not determine AUC, calibration-slope, O/E or clinical
utility precision; those require their own targets and predictive-distribution
assumptions. There is no universal event-count certificate for all metrics.

## Recruitment accounting and the censoring stop

Expected jointly evaluable count is N_candidate*q_joint*q_known only when
q_known is explicitly CONDITIONAL on joint availability. Thus
N_candidate approximately ceil(n_required/(q_joint*q_known)).
This is expectation-based accounting under the stated selection/MCAR assumptions,
not assurance of reaching n, correction for informative missingness, or validation
of the source population. Do not multiply separate marginal assay availability
rates as though independent. Do not reduce the outcome rate by q a second time:
pi must refer to the actual evaluable target domain.

Actual random availability means the above ceiling does not guarantee at least
n; an optional future recruitment-assurance calculation could use an explicitly
assumed Binomial(N,q_eval), not pretend the expectation is an assurance bound.
Uncertain/informative availability or endpoint linkage changes the estimand and
needs a separate missingness/selection strategy, not just more recruitment.

**If c>0 before8years, set `IPCW_method_extension_required` and NO current
uncensored-validation eligibility/adequacy flag.** Do not treat unknown survivors
as noncases, drop them silently, or claim that dividing N by(1-c) validates the
full-horizon estimand. An accounting-only complete-case number may still be shown
but must be distinctly labelled and cannot enable validation. Unknown early
underlying cause also needs its own ascertainment/missing-cause plan.

An IPCW extension would change the loss and its variance. With a positive lower
censor-survival bound g_min, weighted loss differences can be as large as
1/g_min rather than1; estimation of G, weights, ties and clustering contributes
uncertainty. Supplying c alone does not specify G(t), censor/outcome dependence
or the paired IPCW loss variance. Do not reuse the bounded[-1,1] variance/SE as
though that extension were already implemented and qualified.

## Event counts: exact independent-binomial quantiles are not DEFF counts

Under an EXPLICIT independent Bernoulli event assumption:
K~Binomial(n_evaluable,pi), E(K)=n*pi;
report inverse-CDF quantiles at.025/.975 (discrete predictive range), optionally
P(K<an explicitly supplied operational event floor). These are predictive count
quantiles given the assumed pi, NOT a Clopper–Pearson CI for unknown incidence,
not simultaneous limits for both causes, and not power/adequacy evidence.

Do not replace n by n/DEFF in Binomial while still calling results raw events.
Expected raw events remain n*pi. DEFF_D for loss precision need not equal
DEFF_Y for event counts. If clustered event dependence is assumed, exact
binomial quantiles no longer apply: either flag them as inapplicable, simulate
a specified cluster-event process, or separately label the rough normal count
variance DEFF_Y*n*pi*(1-pi) and its limitations/boundaries. There is no exact
DEFF-adjusted binomial distribution without a generative dependence model.
The former50-event feasibility floor is not a sample-size power promise.

## Optional simulation: numerical checks only under declared toy laws

Minimal implementation can simulate bounded D directly on explicit finite
support in[-1,1] with supplied probabilities. A joint scenario may specify
Y~Bernoulli(pi), D|Y=0 and D|Y=1 finite-support laws. Exact moments are:
mu=(1-pi)*mu0+pi*mu1;
v=(1-pi)*v0+pi*v1+pi*(1-pi)*(mu1-mu0)^2.
Each conditional variance must satisfy v_y<=1-mu_y^2. These are stipulated loss
laws, not generated evidence that an implementable assay/model achieves them.

For an optional bounded clustering check with equal cluster size m: for each
cluster, with probability rho draw one common (Y,D) pair for all members,
otherwise draw iid pairs. Marginals stay unchanged and within-cluster Corr(D)
(and Corr(Y), when variable)=rho, giving DEFF=1+(m-1)*rho. Require integer full
clusters and0<=rho<=1. This validates that specific positive-correlation model;
it does not validate arbitrary survey DEFF, algorithm training variability,
clinical performance or true confidence-interval coverage in a real cohort.

Use a fixed seed/replicate count before running. Compare simulated mean/variance
of Dbar with the analytical moments and report Monte Carlo SEs. If assessing
nominal CI coverage, report its Monte Carlo uncertainty; do not tune the seed,
scenario or tolerance to manufacture95%. Simulated direct-loss intervals are
not bootstrap intervals for the33B trained algorithm.

## Implementation and testing boundaries after root protocol

No repo code until the method/scenarios are committed. Proposed pure modules:
scenario validation; analytical precision and recruitment accounting; independent
binomial predictive summaries; optional finite-support simulation; aggregate
report with methods/assumptions/limitations and source/version hashes. No raw
cohort/model-result parser, new source-data download, fit, split search or
clinical recommendation generator.

Meaningful tests:
- exact algebra/sign of paired losses; range/moment feasibility, invalid bools,
  NaN/Inf, probabilities, SD-vs-SE anchoring and q_joint conditional semantics;
- normal inverse/ceiling boundaries; h halves ->n approximately quadruples,
  DEFF doubles ->n approximately doubles, effect sign alone cannot change
  fixed-variance precision N; no hidden observed-power output;
- supplied finite-support moments match exact weighted sums and total variance;
  bounded rare-event counterexample invalidates false pi-only variance limits;
- exact binomial CDF/PPF inversion/boundaries/pi0/1; no DEFF raw-trial substitution;
- recruitment multiplier is expectation-only, q0 impossible, q1 no inflation,
  correlated assay marginals do not masquerade as joint availability;
- ANY c>0 yields IPCW-extension-required and never an uncensored eligible result;
  assay/source/artifact eligibility remains separately blocked when unqualified;
- algorithm_replication always reports unqualified training-level precision;
  locked-source transport eligibility cannot pass without a separately frozen
  pre-external-outcome artifact/provenance and measurement mapping;
- deterministic seed/output aggregate hashes; simulation-moment discrepancies
  judged against predeclared Monte Carlo SE bounds, not brittle lucky draws;
- no reads of33B result estimates, participant files or restricted/network data,
  no fitting/deployment, no automatic sufficient-power/clinical-readiness claims.

## Primary methodological anchors checked for this design

Riley et al.,2021, Statistics in Medicine, DOI10.1002/sim.9025:
https://pubmed.ncbi.nlm.nih.gov/34031906/ . Their external-validation planning
framework targets prespecified precision of distinct performance measures and
requires model/population assumptions; it does not supply a generic events rule
or establish the paired-loss formula above. The formula here is an elementary
variance calculation for a mean of paired loss differences with an explicitly
assumed DEFF, not attributed to that article as a validated extension.

SciPy official binom documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom.html .
PPF is the inverse CDF of the independent n-trial binomial. Pin the actual
implementation version later; reading current documentation is not a dependency
upgrade. Only public methodological documentation was browsed; no cohort records,
new outcome estimates, account actions or access agreements were requested.
