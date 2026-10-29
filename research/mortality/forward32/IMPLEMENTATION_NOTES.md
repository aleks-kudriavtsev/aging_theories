# Pre-result implementation clarifications

These notes clarify execution of plan377778c before any new temporal metrics.
They change no split, target, predictor, weight, horizon, penalty, contrast,
bootstrap mode/count/seed or source interpretation.

The source adapter exports the three-fold mapping derived from eligible early
clusters. The numerical fit helper computes the same hash-order mapping from
common early clusters. To prove equivalence rather than silently change other
folds when a whole cluster is lost, the feasibility gate additionally requires
common and eligible masked-PSU inventories to be identical in EACH period.
Failure produces blocked_feasibility and no fits. With equal cluster sets,
hash-sort/rank modulo3 mapping is identical; retaining zero-weight rows in
resampling preserves that mapping. Synthetic tests compare the exported mapping
with an independent SHA256 implementation and check lost-cluster rejection.

Calendar separation refers to biological BASELINE cohorts. The eight-year
early labels become known in2007–2010, after the late2003–2004 baselines, and
the common public mortality release is2019. Thus this is retrospective baseline-
calendar stress testing, not a historical as-of2003 model or prospective
deployment validation. No late participant's outcome enters fitting or tuning;
the stated label-availability limitation still matters to interpretation.

For descriptive feasibility, later labels are used only to report support and
apply the fixed go/block threshold. Passing late-label changes cannot change
the earlier fitted parameters. Later calibration regression is diagnostic only
and never changes evaluated predictions. Both uncertainty modes remain separate;
joint_refit is the primary uncertainty and fixed_training excludes training
uncertainty. Failures remain local to arms/metrics with no redraw.
