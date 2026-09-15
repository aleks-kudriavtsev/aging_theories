# Iteration23: NfL robustness and joint protein source audit

Date:2026-09-15. Base cb19a2a2c7e5e1cc42e008595453ba2a51acf24c.

This is an explicitly post-hoc robustness study. The prior NfL result (1270 participants,37 deaths within4years; unfavorable incremental Brier) is already known. No new independent validation is claimed. The previous protocol, model, outputs and public calculator remain unchanged.

## Frozen scope
Reuse the exact official2013-2014 CDC files from run34994970265/artifact10406939200 and verify their original manifest. Reconstruct the old complete domain and old outcome-blind whole-PSU fold allocation. Baseline expanded13 remains fixed (SHA25687a8623b9fd5cb793ce2da5311a17a56b7d94251da5df335de3d406e437a55da). Reproduce the previous OOF estimates before new analysis.

## Added uncertainty analysis
1000 paired rescaled PSU bootstrap draws, RNG seed20260923. Draw in the broader eligible sample, restrict to the same complete domain, then REFIT calibration-only and calibration+log2NfL within each original training fold using replicate survey weights. No participant or copies of the same PSU cross the fixed fold boundary. Same ridge1 and original weight-normalization algorithm. Score both OOF models with the same replicate weights. Report percentile intervals, undefined replicate/fold counts and paired comparison with the fixed-OOF approach using the identical draws. This adds uncertainty from estimating the small NfL update conditional on the frozen baseline, original fold allocation and selected population; it does not cover baseline development uncertainty, model-selection history or external transport.

## Partition sensitivity
Repeat the original OOF algorithm for20 additional deterministic outcome-blind PSU partitions (salt joint23:split:0 through19, SHA256 sorted then round-robin5). Retain all, report signed differences and range; do not select a favorable partition or call repetitions independent validations. No tuning, alternative biomarker forms, cutoff search or subgroup-based model selection.

## Joint-data acquisition
Resolve exact released variables and measurement forms for GDF15,IL6,TNF,GFAP,NfL,NT-proBNP and cystatinC against public HRS descriptions/codebooks. Distinguish available metadata, authorized raw records and actual same-person/same-baseline intersection. Never substitute sTNFR1 for TNF or invent absent UACR/HbA1c. Protected data and signatures are not obtained or published without specific authorization. Public source files may be fingerprinted and audited; no mortality training from sample data lacking mortality or the full baseline panel.

## Release
Publish code, aggregate metrics, source metadata and checks only. No participant IDs or records. Clinical readiness remains false. Any deviation from this plan is documented without replacing the plan.
