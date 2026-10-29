"""Scenario-only conditional paired-loss precision; no data/model imports.

Arithmetic and simulations condition on stipulated loss/prediction laws. They
do not qualify clinical validation, nested-algorithm precision or sample power.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.stats import binom, norm

MAX_COUNT = 2**53 - 1
_REQUIRED = frozenset({'validation_mode', 'mu_d', 'alpha', 'half_width', 'n_evaluable',
    'deff_d', 'deff_basis', 'q_joint', 'q_known_conditional', 'censoring_fraction',
    'raw_evaluable_event_probability', 'event_count_model', 'deff_event'})
_OPTIONAL = frozenset({'id', 'sd_d', 'se_reference', 'n_reference', 'deff_reference',
                       'technical_event_floor'})


def _number(value, name):
    if type(value) not in (int, float):
        raise ValueError('Finite numeric value required: ' + name)
    try:
        valid = math.isfinite(value)
    except (OverflowError, ValueError):
        valid = False
    if not valid:
        raise ValueError('Finite numeric value required: ' + name)
    return float(value)


def _count(value, name, minimum=1):
    if type(value) is not int or not minimum <= value <= MAX_COUNT:
        raise ValueError('Representable integer count required: ' + name)
    return value


def _probability(value, name, open_interval=False):
    value = _number(value, name)
    if not (0 < value < 1 if open_interval else 0 <= value <= 1):
        raise ValueError('Probability outside allowed range: ' + name)
    return value


def _positive(value, name):
    value = _number(value, name)
    if value <= 0:
        raise ValueError('Positive value required: ' + name)
    return value


def _ceil_count(value, name):
    if not math.isfinite(value) or value > MAX_COUNT:
        raise ValueError('Planning count exceeds numerical representation: ' + name)
    return max(1, math.ceil(value))


def plan(scenario):
    """Return conditional arithmetic only, from explicit scenario assumptions.

    Exactly one variance source is accepted: sd_d OR the complete anchored
    se_reference/n_reference/deff_reference triple. No observed-result parser,
    default effect size, independent marginal assay product or power output.
    """
    if (not isinstance(scenario, dict) or not _REQUIRED <= set(scenario)
            or set(scenario) - (_REQUIRED | _OPTIONAL)):
        raise ValueError('Exact explicit precision-scenario fields required')
    if 'id' in scenario and (not isinstance(scenario['id'], str) or not scenario['id']):
        raise ValueError('Nonempty scenario identifier required')
    mode = scenario['validation_mode']
    if mode not in ('algorithm_replication', 'fixed_scorer'):
        raise ValueError('Unknown validation conditioning')
    if scenario['deff_basis'] != 'total_relative_to_iid_unweighted_n':
        raise ValueError('One total DEFF required; no implicit Kish/design double counting')
    mu = _number(scenario['mu_d'], 'mu_d')
    if not -1 <= mu <= 1:
        raise ValueError('Paired Brier-loss mean outside[-1,1]')
    anchored = {'se_reference', 'n_reference', 'deff_reference'}
    if 'sd_d' in scenario:
        if anchored & set(scenario):
            raise ValueError('Choose paired SD or anchored SE, not both')
        sd = _number(scenario['sd_d'], 'sd_d')
        if sd < 0:
            raise ValueError('Paired SD must be nonnegative')
        variance = sd * sd
        positive_scale = sd > 0
        variance_source = 'specified_individual_paired_loss_SD'
    else:
        if not anchored <= set(scenario):
            raise ValueError('Reference SE needs both reference N and reference DEFF')
        se_ref = _number(scenario['se_reference'], 'se_reference')
        if se_ref < 0:
            raise ValueError('Reference SE must be nonnegative')
        n_ref = _count(scenario['n_reference'], 'n_reference')
        deff_ref = _positive(scenario['deff_reference'], 'deff_reference')
        variance = se_ref * se_ref * n_ref / deff_ref
        positive_scale = se_ref > 0
        variance_source = 'anchored_SE_inverted_under_conditional_loss_variance_assumption'
    bound = (1 - abs(mu)) * (1 + abs(mu))
    if positive_scale and variance == 0:
        raise ValueError('Positive variance underflows numerical representation')
    feasible_roundoff = math.isclose(variance, bound,
        rel_tol=8 * np.finfo(float).eps, abs_tol=0.)
    if (not math.isfinite(variance)
            or (variance > bound and not feasible_roundoff)):
        raise ValueError('Mean/variance incompatible with D in[-1,1]')
    clipped = variance > bound
    variance = min(variance, bound)
    alpha = _probability(scenario['alpha'], 'alpha', open_interval=True)
    half_width = _positive(scenario['half_width'], 'half_width')
    n = _count(scenario['n_evaluable'], 'n_evaluable')
    deff = _positive(scenario['deff_d'], 'deff_d')
    z = float(norm.isf(alpha / 2))
    mean_variance = deff * variance / n
    se = math.sqrt(mean_variance)
    ratio = z * math.sqrt(deff * variance) / half_width
    required_n = _ceil_count(ratio * ratio, 'conditional_n_for_half_width')
    q_joint = _probability(scenario['q_joint'], 'q_joint')
    q_known = _probability(scenario['q_known_conditional'], 'q_known_conditional')
    censoring = _probability(scenario['censoring_fraction'], 'censoring_fraction')
    yield_probability = q_joint * q_known
    recruitment = (_ceil_count(required_n / yield_probability, 'candidate_recruitment_expectation')
                   if yield_probability > 0 else None)
    pi = _probability(scenario['raw_evaluable_event_probability'], 'raw_evaluable_event_probability')
    event_deff = _positive(scenario['deff_event'], 'deff_event')
    event_model = scenario['event_count_model']
    if event_model not in ('independent_binomial', 'dependent_unspecified'):
        raise ValueError('Explicit independent or unspecified-dependent event model required')
    if event_model == 'independent_binomial' and event_deff != 1:
        raise ValueError('Exact binomial counts require event independence/DEFF_event1')
    floor = scenario.get('technical_event_floor')
    if floor is not None:
        floor = _count(floor, 'technical_event_floor', minimum=0)
    iid = event_model == 'independent_binomial'
    predictive_range = ([int(binom.ppf(alpha / 2, n, pi)), int(binom.ppf(1 - alpha / 2, n, pi))]
                        if iid else None)
    probability_below = (float(binom.cdf(floor - 1, n, pi)) if iid and floor is not None else None)
    return dict(schema='replication34-conditional-precision/1', scenario_id=scenario.get('id'),
        validation_mode=mode, mean_d=mu, paired_loss_variance=variance,
        paired_loss_SD=math.sqrt(variance), variance_source=variance_source,
        moment_variance_roundoff_clipped=clipped, alpha=alpha, normal_quantile=z,
        n_evaluable=n, deff_d=deff, deff_basis=scenario['deff_basis'],
        conditional_mean_variance=mean_variance, conditional_SE=se,
        conditional_normal_half_width=z * se, requested_half_width=half_width,
        conditional_n_for_half_width=required_n, mean_effect_used_to_choose_n=False,
        degenerate_zero_loss_variance=variance == 0,
        q_joint=q_joint, q_known_conditional=q_known,
        joint_evaluable_yield_probability=yield_probability,
        candidate_recruitment_expectation=recruitment,
        recruitment_target_reachable_under_assumed_availability=yield_probability > 0,
        recruitment_is_expectation_not_assurance=True,
        recruitment_MCAR_accounting_assumption_not_selection_correction=True,
        censoring_fraction=censoring, IPCW_method_extension_required=censoring > 0,
        bounded_uncensored_loss_formula_applies_to_full_CIF=censoring == 0,
        full_CIF_validation_eligible=False,
        conditioning='conditional_fixed_prediction_loss_proxy_not_full_refitted_algorithm',
        algorithm_precision_status='adequacy_unqualified',
        algorithm_training_uncertainty_included=False,
        frozen_source_artifact_required_before_external_outcomes=(mode == 'fixed_scorer'),
        artifact_or_assay_qualification_performed=False,
        normal_interval_is_approximation_not_coverage_guarantee=True,
        clinical_adequacy_or_power_established=False,
        calibration_or_AUC_precision_established=False,
        event_count=dict(model=event_model, raw_evaluable_probability=pi,
            n_raw_trials=n, expected_raw_events=n * pi,
            variance_proxy=event_deff * n * pi * (1 - pi), deff_event=event_deff,
            predictive_range=predictive_range, predictive_range_exact_distribution_assumption=iid,
            predictive_range_is_not_incidence_CI=True,
            technical_floor=floor, probability_below_technical_floor=probability_below,
            technical_floor_establishes_power_or_adequacy=False,
            n_over_DEFF_used_as_binomial_trials=False))


def _atoms(atoms):
    if not isinstance(atoms, list) or not atoms:
        raise ValueError('Nonempty finite prediction-atom list required')
    checked = []
    for row in atoms:
        if not isinstance(row, dict) or set(row) != {'q', 'r', 'p_new', 'p_ref'}:
            raise ValueError('Exact prediction-atom fields required')
        checked.append([_probability(row[name], name) for name in ('q', 'r', 'p_new', 'p_ref')])
    values = np.asarray(checked, float)
    total = float(values[:, 0].sum())
    if not math.isclose(total, 1., rel_tol=0., abs_tol=1e-12):
        raise ValueError('Prediction-atom probabilities must sum to1')
    values[:, 0] /= total
    return values, total


def finite_prediction_atoms_moments(atoms):
    """Exact moments: draw X atom, fix both predictions, THEN draw Y|X."""
    values, total = _atoms(atoms)
    q, r, p_new, p_ref = values.T
    difference = p_new - p_ref
    d0, d1 = difference * (p_new + p_ref), difference * (p_new + p_ref - 2)
    conditional_mean = (1 - r) * d0 + r * d1
    mean = float(q @ conditional_mean)
    second = float(q @ ((1 - r) * d0**2 + r * d1**2))
    within = float(q @ (4 * difference**2 * r * (1 - r)))
    between = float(q @ (conditional_mean - mean)**2)
    variance = max(0., within + between)
    return dict(schema='replication34-prediction-atoms-moments/1', atoms_n=len(atoms),
        input_probability_sum=total, normalization_roundoff_applied=total != 1.,
        raw_event_probability=float(q @ r), mean_d=mean, second_moment_d=second,
        variance_d=variance, SD_d=math.sqrt(variance),
        mean_conditional_Y_variance=within, variance_of_conditional_means=between,
        variance_decomposition_error=abs(second - mean * mean - variance),
        predictions_fixed_before_Y=True, fitted_model_or_clinical_validation=False)


def _coverage(count, total):
    value = count / total
    z = float(norm.isf(.025))
    denominator = 1 + z*z/total
    center = (value + z*z/(2*total)) / denominator
    radius = z * math.sqrt(value*(1-value)/total + z*z/(4*total*total)) / denominator
    return dict(covered=count, replicates=total, proportion=value,
        binomial_MCSE=math.sqrt(value * (1 - value) / total),
        MC_95pct_Wilson_interval=[max(0., center-radius), min(1., center+radius)])


def simulate_case(atoms, n_full_clusters, cluster_size, rho, replicates, seed, alpha=.05):
    """Toy-law simulation only; independent full clusters, shared-pair mixture.

    With probability rho one (X,Y) pair is shared within a cluster; otherwise
    members are iid. Marginals are preserved and DEFF=1+(m-1)*rho exactly.
    No fitted prediction algorithm or external cohort is simulated/validated.
    """
    j = _count(n_full_clusters, 'n_full_clusters')
    m = _count(cluster_size, 'cluster_size')
    n = _count(j * m, 'total_simulated_participants')
    repetitions = _count(replicates, 'replicates', minimum=2)
    if type(seed) is not int or seed < 0:
        raise ValueError('Explicit nonnegative integer seed required')
    rho = _probability(rho, 'rho')
    alpha = _probability(alpha, 'alpha', open_interval=True)
    values, _ = _atoms(atoms)
    exact = finite_prediction_atoms_moments(atoms)
    q, risk, p_new, p_ref = values.T
    deff = 1 + (m - 1) * rho
    theory_variance = deff * exact['variance_d'] / n
    theory_se = math.sqrt(theory_variance)
    z = float(norm.isf(alpha / 2))
    rng = np.random.default_rng(seed)
    means = np.empty(repetitions); counts = np.empty(repetitions)
    known_coverage = estimated_coverage = 0
    for rep in range(repetitions):
        indices = rng.choice(len(q), size=(j, m), p=q)
        y = rng.random((j, m)) < risk[indices]
        shared = rng.random(j) < rho
        if shared.any():
            shared_indices = rng.choice(len(q), size=int(shared.sum()), p=q)
            shared_y = rng.random(len(shared_indices)) < risk[shared_indices]
            indices[shared] = shared_indices[:, None]
            y[shared] = shared_y[:, None]
        a, b = p_new[indices], p_ref[indices]
        d = (a - b) * (a + b - 2*y)
        cluster_means = d.mean(axis=1)
        means[rep], counts[rep] = float(cluster_means.mean()), int(y.sum())
        deviation = abs(means[rep] - exact['mean_d'])
        if theory_variance == 0:
            #Only a declared degenerate law permits roundoff handling. Scaling
            #this check to its loss magnitude cannot widen a positive-variance
            #CI or turn tiny but genuine sampling variation into coverage.
            loss_scale = max(abs(exact['mean_d']), float(np.max(np.abs(d))))
            roundoff_bound = 16 * np.finfo(float).eps * loss_scale
            if (deviation > roundoff_bound
                    or np.max(np.abs(d - exact['mean_d'])) > roundoff_bound):
                raise ArithmeticError('Zero-variance toy law produced nonconstant losses')
            known_coverage += 1
            estimated_coverage += int(j > 1)
        else:
            known_coverage += int(deviation <= z * theory_se)
            if j > 1:
                estimated_se = math.sqrt(float(cluster_means.var(ddof=1)) / j)
                estimated_coverage += int(deviation <= z * estimated_se)
    observed_variance = float(means.var(ddof=1))
    centered_fourth = float(np.mean((means - means.mean())**4))
    variance_mcse = math.sqrt(max(0., (centered_fourth
        - ((repetitions-3)/(repetitions-1))*observed_variance**2) / repetitions))
    pi = exact['raw_event_probability']
    return dict(schema='replication34-toy-simulation/1', seed=seed, replicates=repetitions,
        n_full_clusters=j, cluster_size=m, n_participants=n, rho=rho, deff=deff, alpha=alpha,
        exact_moments=exact, exact_mean_d=exact['mean_d'], exact_variance_mean_d=theory_variance,
        exact_normal_half_width=z*theory_se,
        MC_mean_d=float(means.mean()), MC_variance_mean_d=observed_variance,
        MC_centered_fourth_moment_dbar=centered_fourth,
        MCSE_mean_d=theory_se / math.sqrt(repetitions),
        estimated_MCSE_mean_d=math.sqrt(observed_variance / repetitions),
        approximate_fourth_moment_MCSE_variance=variance_mcse,
        known_toy_variance_normal_CI_coverage=_coverage(known_coverage, repetitions),
        estimated_cluster_mean_SE_normal_CI_coverage=(_coverage(estimated_coverage, repetitions) if j > 1 else None),
        positive_variance_CI_coverage_uses_no_added_tolerance=True,
        exact_zero_variance_handling='assert_constant_losses_within16eps_of_loss_scale_then_cover_degenerate_mean',
        exact_expected_event_count=n*pi,
        exact_event_count_variance_in_this_toy=deff*n*pi*(1-pi),
        MC_mean_event_count=float(counts.mean()), MC_variance_event_count=float(counts.var(ddof=1)),
        MCSE_mean_event_count=math.sqrt(deff*n*pi*(1-pi)/repetitions),
        event_and_loss_DEFF_equal_only_under_this_specified_shared_pair_law=True,
        simulation_scope='stipulated_prediction_atoms_and_cluster_mixture_only',
        predictions_fixed_before_Y=True, refitted_algorithm_precision_validated=False,
        external_cohort_or_clinical_validation_performed=False,
        participant_simulation_arrays_exported=False)
