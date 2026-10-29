# Independent Panel25 audit, 2026-10-29

Read-only audit; no ZIP download, no empirical recomputation, no repository files changed.

## Verdict

No material error found in the declared exploratory comparison. The retained conclusion is warranted: linear log2 serum urate did not demonstrate incremental four-year mortality prediction over calibration-only in this reused NHANES 2013–2014 complete-case MEC domain. It does not establish absence of biological or prognostic value of urate generally. The study cannot be described as new independent validation or an overall model improvement.

## Independent aggregate and fingerprint checks

Obtained original job123474532324/run40881679826 and already-completed repeat job123475996644/run40882094811 through authorized github_fetch_workflow_job_logs. Parsed the JSON stdout directly, without downloading artifacts.

- Local RESULTS_SUMMARY.json is structurally identical to the original job JSON, every field.
- Both jobs: n2867, events4=96, eligible_n3218, complete_baseline_n2868; source survey10175 and age40–75 count3308. Flows add correctly:3218−350−1=2867.
- All reported AUCs and all bootstrap interval limits/counts match exactly between the two existing jobs. No failed/undefined bootstrap refits in either job.
- Original ΔBrier0.000004504522136319217; repeat0.000004504523340727318. Absolute difference1.2044081010298413e−12, negligible optimizer/platform variation. The runs are numerically reproducible, not byte-identical result JSONs.
- Local analyze.py SHA256 c0188b3eb4255bc01a1c24c0dd60d9a7a13c6dfc1d116a174584a7b089c4bd98, plan1cd59357168dde9113bff46d311dde535201373303568f651f42006776524158 and frozen baseline87a8623b9fd5cb793ce2da5311a17a56b7d94251da5df335de3d406e437a55da match both jobs exactly.
- Mortality source SHA588acea1d4cdf47a14cc9a35bd56e372739b8788d63ca0af98bf954c70c29905 matches both jobs; both independent readers report81400 exact cell matches.
- Source manifest SHA differs: original0de17151efc7cf0baa3292d8bddb3ef47f5d456d15a549efd8b9c0d84655b0a9, repeat11ad6bb7a6b3d2d01ad11a24ac36e0a22f331b617892af845563c8fbce67538b. retrieve.py embeds retrieval date_utc, so manifest SHA is inherently variable. The stdout does not expose individual XPORT fingerprints, therefore equality of all downloaded scientific source bytes remains unverified. Do not infer it from aggregate metric agreement.
- Local summary byte SHA c9bf4bef5c8827f5357ce4e66139f6bc84a09ed81a90e23b5e5aaac8dad82e80. ZIP digest has not been independently verified; no attempt made to bypass prior403.

## Scientific design

1. Likelihood is appropriate for an all-cause multiplicative update of the fixed summed cause-specific hazard: weighted D*(a+b*x)−H0(min(T,4))*exp(a+b*x), plus ridge1 on b. Omitted baseline log hazard term is constant in update coefficients. Four years lies wholly in baseline interval0–5; no interval extrapolation or Cox fit occurs.
2. Frozen predictor transformations and coefficients are never fit on Panel25. Training-only centering of urate and calibration fitting occurs within whole-PSU folds. assign_folds reads only cluster identities; no outcome-driven assignment found. Complete-case masks are shared across both compared updates. The existing split test checks that a fold's outcomes cannot alter its own predictions.
3. WTMEC2YR is a defensible weight for the broader examined cohort/routine laboratory panel. WTSSNH2Y would wrongly restrict/weight the NfL subsample; Panel25 correctly avoids it. Inference still targets mortality-eligible complete profiles, with no nonresponse/complete-case bias correction claimed.
4. Alive participants with t_exam<48 cause a hard stop before metrics. Deaths at month0 receive0.5months for likelihood exposure; death status at four years remains unchanged. For this release, binary weighted Brier/AUC does not need IPCW after that gate.
5. Bootstrap is paired at PSU within strata, generated on eligible MEC before complete-case restriction. Each update is refit using replicate weights with original fold boundaries. Counts of undefined statistics are explicitly retained. Intervals remain conditional on a previously developed baseline, the chosen panel, functional form and partition; they do not cover baseline training/history/external transport.

## Concrete limitations and improvements for next iteration

Priority1 — Provenance coverage: Panel25 records only analyze.py code hash although fitting/scoring/readers come from imported nfl22/analyze.py, joint23/robustness.py, competing_risks09.py, transport_no_crp12.py, validate_temporal08.py and xpt_checked.py. Record dependency-file hashes and git commit plus a stable source-scientific-file fingerprint map in future result outputs. A code hash alone is insufficient to replay the scientific computation.

Priority1 — Auditability of real fits: urate_crossfit discards fold fit rows. Persist only aggregate fold train/test n/events, parameter/centering/convergence diagnostics, strata/PSU counts, weighted effective sample size, event-weight support and minimum survivor followup. No individual rows/IDs/predictions needed. These allow real-data leakage/support/convergence review beyond merely successful runner exit.

Priority2 — Numerical fail-closed guard inherited from NfL22: fit_update checks finite coefficients and gradient>threshold but does not explicitly reject a nonfinite gradient or objective; nan>threshold is false. Existing runs contain no demonstrated failure, so this is a defensive next-iteration improvement, not grounds to reject observed estimates. Validate finite objective/gradient and preserve diagnostics without rewriting the frozen implementation/results.

Priority2 — Archive-derived baseline source: load_model verifies exact frozen JSON, but the original Integrated20/model_bundle20.json is absent here. Its original extraction equality/provenance cannot be reverified from this checkout, and Panel25 does not claim such a fresh parent verification. Preserve original authoritative provenance documentation and avoid attributing its training history to a newly verified parent.

Priority2 — New cohorts: this audit does not establish that any later cycle is untouched. Inspect all remote branches, scientific scripts/reports and prior supplied-archive descriptions before calling a cohort independent. If the true history is unknown, classify new calculation as validation on a previously investigated population. Do not choose alternative urate transformations based on current negative result.

No blocker prevents a separately planned public-data iteration. Only the previously denied local ZIP byte comparison and full scientific-source-map comparison remain unavailable; already-completed job stdout suffices for the limited independent aggregate verification above.
