# Stage34 paired-Brier precision method checks (read-only)

Pure mathematical review. No repository edits, participant/new outcome reads,
source GET, model fits or empirical effects. In particular, no33B effect or CI
is used as a default or as post-hoc power. No clinical threshold is proposed.

## 1. Paired difference, sign and variance

Let a=p_new, b=p_ref, Y in{0,1}, with both probabilities in[0,1]. Define
D=(Y-a)^2-(Y-b)^2, so negative meanD means lower Brier for the new model.
With delta=a-b and s=a+b,

    D=delta*(s-2Y).

Each squared loss is in[0,1], hence D in[-1,1]. If mu=E[D],

    0<=v=Var(D)<=1-mu^2<=1.

This is sharper than a variance bound1 when a nonzero mu is specified. Any
scenario with |mu|>1 or v<0 or v>1-mu^2 is invalid. SD, variance and SE are not
interchangeable inputs: SD=sqrt(v), SE approximatelysqrt(DEFF_D*v/n_eval).

For fixed covariates/prediction pair X, let r(X)=P(Y=1|X). Then

    E[D|X]=delta*(s-2r),
    Var(D|X)=4*delta^2*r*(1-r),
    Var(D)=E[4*delta^2*r*(1-r)] + Var(delta*(s-2r)).

The second term is essential. Var(D)<=4*pi*(1-pi) is NOT a general bound from
rare-event prevalence pi. Counterexample: equally likely prediction pairs
(a,b)=(1,0) and(0,1), independent of Y~Bernoulli(pi). D is symmetric{-1,+1},
mean0 and variance1 for ANY pi. Event prevalence alone cannot identify paired
score variance, error correlation or the precision of incremental Brier.

If a scenario also supplies mean Briers A=E[(Y-a)^2], B=E[(Y-b)^2], then
mu=A-B and another valid bound is

    v <= min(A+B, 2-A-B) - mu^2.

Proof: for losses Lnew,Lref in[0,1], D^2<=|D|<=Lnew+Lref and
D^2<=|D|<=2-Lnew-Lref. Such marginal losses must themselves be valid;
this algebra is not evidence that a particular risk model achieves them.

If |a-b|<=epsilon<=1 is an explicit scenario restriction, then
|D|<=2epsilon-epsilon^2, a tighter geometry bound than2epsilon. The conditional
noise bound uses4delta^2*r*(1-r), but total variation still needs the mean term.

## 2. Given-Y consistency and honest bounded scenario simulation

For class-conditional means mu0=E[D|Y=0], mu1=E[D|Y=1] and variances v0,v1,

    mu=(1-pi)*mu0+pi*mu1,
    v=(1-pi)*v0+pi*v1 + pi*(1-pi)*(mu1-mu0)^2,
    -1<=mu_y<=1, 0<=v_y<=1-mu_y^2.

The between-class term gives a LOWER bound on totalv once both conditional
means and pi are fixed. Do not allow inconsistent independent choices of
marginalmu/v and conditionalmu_y/v_y. At pi=0/1, the absent class is unidentified;
its conditional moments should not be treated as empirically known.

An exact bounded-D moment simulator, for an explicitly hypothetical D-only
scenario, is:
- for |mu_y|<1 set alpha_y=v_y/(1-mu_y^2);
- with probability1-alpha_y emit D=mu_y;
- with probabilityalpha_y emit+1 with probability(1+mu_y)/2, otherwise-1;
- for |mu_y|=1 requirev_y=0 and emit that constant.
Draw Y~Bernoulli(pi), then its specified conditionalD distribution. This realizes
the stated moments exactly, never emits an impossible D outside[-1,1], and does
not invent empirical effect estimates. It studies a bounded score-difference
sampling distribution; it is NOT a simulation/validation of a fitted predictive
algorithm or clinical risk model.

A Gaussian draw followed by clipping is not an equivalent implementation: it
changes the mean/variance and may create endpoint point masses. An unbounded
Gaussian D draw is incompatible with the score support.

For a prediction-level simulator without outcome-derived predictions, prefer
finite X atoms j: fixed probabilities q_j, true risks r_j, and prediction pairs
a_j,b_j all in[0,1]. Draw j~Categorical(q), set predictions(a_j,b_j), THEN draw
Y~Bernoulli(r_j). Analytical audit:

    pi=sum(q_j*r_j),
    m_j=(a_j-b_j)*(a_j+b_j-2r_j),
    mu=sum(q_j*m_j),
    v=sum(q_j*(4*(a_j-b_j)^2*r_j*(1-r_j)+m_j^2))-mu^2.

GivenY, P(j|Y=1)=q_j*r_j/pi and P(j|Y=0)=q_j*(1-r_j)/(1-pi), so conditional
moments can be verified as well. For a SINGLE fixed prediction pair, D0=delta*s
and D1=delta*(s-2) are coupled: s in[|delta|,2-|delta|]. They cannot be chosen
independently. In particular, a>b harms the Y=0 squared loss and helps Y=1.
Do not back-solve predictions using actual observedY then present the resulting
numbers as deployable-model simulation; that would encode the label.

## 3. Approximate precision, not a significance/power guarantee

For fixed prespecified models and a paired evaluation mean under a declared
sampling design,

    Var(meanD) approximately DEFF_D*v/n_eval,
    halfwidth approximately z_(1-alpha/2)*sqrt(DEFF_D*v/n_eval),
    n_eval approximately ceil(z^2*DEFF_D*v/h^2).

Here n_eval is actual evaluable paired rows and DEFF_D is the design effect for
THIS paired score difference relative to the declared IID baseline. Positive
h, finitevalidv, positiveDEFF and an explicit intervallevel are required. The
normal approximation is not exact finite-sample coverage, especially with few
independent PSUs, rare classes, skewed differences or model selection. If a
survey t approximation is intended, its degrees of freedom/method must be
specified rather than mechanically applying z or an unreviewed PSU-minus-strata
rule to a complex refit procedure.

A totalDEFF may include unequal weights and cluster correlation. Do not first
replace n by Kish n_eff=(sumw)^2/sumw^2 and ALSO multiply by a DEFF that already
includes weighting. Under independent equal-varianceD and fixed weights,
Var(weightedmeanD)=v/n_eff; further cluster effects require an additional,
clearly distinct factor. A survey-weighted population incidence is not in
general the raw enrolled-cohort event probability for count planning.

No minimally important clinical effect or decision threshold follows fromh.
Do not infer chance of a favorable result, significance or replication success
from a CIhalfwidth alone. A minimum event/support computational rule, if any,
is a separate input/constraint, not a consequence of the precision equation.

## 4. Recruitment/availability accounting

Let q_joint be the JOINT availability of the required assay/core panel, not the
product of assay marginals. Let q_known=P(endpoint known|panel available). Then
r=q_joint*q_known is the evaluable probability by definition; if q_known is
instead an unconditional marginal, multiplying additionally assumes independence.
Under MCAR / the stated accounting model,

    E[n_evaluable] = N_candidates*r,
    N_candidates approximately ceil(n_eval/r).

This is expected-yield accounting, not a guarantee ofn_eval, representativeness
or an unchanged mu/v after selection. Selection related to outcomes/predictions
can change the estimand and score distribution. q=0 means no finite recruitment
solution, not a numerical divide-by-zero workaround. If only two availability
marginals are known, their intersection is bounded by max(0,q1+q2-1) and
min(q1,q2); independence must not silently become observed overlap.

With independent candidate availability, K~Binomial(N,r). To guarantee
P(K>=n_eval)>=1-alpha_yield under that model, solve its exact binomial tail for
N; N=n_eval/r only fixes the expectation. Assay/endpoint overlap uncertainty
needs separate scenarios or evidence, not a deterministic inflation multiplier.

## 5. Event uncertainty: exact binomial versus design effects

For FIXED evaluable n and IID Bernoulli event probabilitypi,
E~Binomial(n,pi). Exact P(e_min<=E<=n-c_min) is

    BinomCDF(n-c_min;n,pi)-BinomCDF(e_min-1;n,pi).

Events/non-events are complementary, not independent constraints whose
probabilities may be multiplied. If evaluable count is random amongN candidates
and availability is independent of outcome, event/non-event/unavailable counts
are multinomial(N; r*pi,r*(1-pi),1-r). Event marginal is Binomial(N,r*pi), but
event and non-event count tails are dependent. A joint multinomial calculation
or conservative explicitly allocated union bound is needed for joint assurance.

The “exact” label is conditional on these IID/sample-availability assumptions.
With clustered recruitment, unequal sampling or event-dependent missingness,
actual event counts need a design/generative count model. DEFF_D for Brier is
NOT DEFF_Y for event indicators; one cannot substitute n/DEFF into a Binomial
and call that exact. A scalar DEFF specifies a variance ratio, not a full count
law or event-count tail. Expected actual events remain n*pi, not(n/DEFF)*pi.
If a beta-binomial/cluster-risk model is proposed, its sizes/ICC/distribution
must be explicit and separately checked; many count laws share the sameDEFF.
Binomial confidence intervals for weighted survey prevalence are likewise not
“exact” merely because raw event counts exist.

## 6. Censoring or unknown outcomes is a method extension

If some participants are censored beforetau, complete-case deletion is not
automatically an unbiased CIF/Brier assessment. q_known inflation cannot by
itself fix this. The present known-label pairedD formula applies only when the
fixed-horizon binary target is known for every evaluated participant (or an
explicit justified missing-data estimand/assumption has been established).

A censoring>0 path needs a separately specified/validated IPCW or equivalent
survival/competing-risk scoring procedure, independent/conditional censoring and
positivity assumptions, tie handling, and uncertainty from estimating censoring
weights. Known competing deaths are events with targetY=0, not ordinary censoring.
Unknown death cause is a missing-cause problem, not automatically right censoring.

For example, IPCW score contributions use1/G(T-) for observed deaths beforetau
and1/G(tau) for known survivors through tau. Paired weighted contributions can
exceed[-1,1] whenG<1; therefore the raw bounded-D variance/simulator contract
cannot be copied unchanged. Normalization and fittedG can introduce further
variance. A planner supporting only known binary labels should mark censoring>0
as requiring an unimplemented method extension, not declare validity after
dropping records or scaling recruitment by1/(1-c).

## 7. Algorithm replication versus frozen coefficient transport

Frozen-model external evaluation keeps BOTH compared prediction functions,
preprocessing, units/assay mappings and coefficients fixed before evaluation
outcomes. Precision then concerns the evaluation distribution conditional on
those fitted models; no training-uncertainty claim follows. Assay/provenance
qualification is still required; schema compatibility or a unit conversion is
not proof that a historical assay mapping transports.

Algorithm replication refits/tunes a prespecified algorithm in a new cohort and
evaluates with a genuinely separate test set or honest nested cross-fitting.
It is a different estimand from transporting frozen coefficients. Shared fits
make prediction differences dependent; a simple fixed-model n_eval equation
omits training uncertainty unless explicitly conditional on training.

With independent trainingT and evaluationE,

    Var(meanD)=E_T[Var_E(meanD|T)] + Var_T[E_E(meanD|T)]
             approximately DEFF_eval*E_T[v_D(T)]/n_eval + Var_T[mu_D(T)].

Increasing evaluation n alone does not remove the second term. Joint training+
evaluation refit simulation/resampling and explicit n_train/n_eval allocation
are needed for unconditional algorithm-precision planning. Re-estimating an
intercept or assay conversion from evaluation outcomes is a model update, not
unchanged frozen-model validation. Retrofitting after seeing results requires
a new assessment boundary; it must not inherit the old precision/validation
label silently.
