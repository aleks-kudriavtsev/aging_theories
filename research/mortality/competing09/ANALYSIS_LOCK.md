# Competing-risk analysis lock 09

Date: 2026-09-10. Parent: 6b3fe1372c5147aabec6e9e6fc9a87fae8704656. This decision log precedes retrieval/inspection of NHANES 2009-2010 participant outcomes in this iteration. It is not an independent registry entry. NHANES2005-2008 has already been evaluated for all-cause mortality and is now explicitly development data, not a new holdout.

## Question and endpoint partition
Test whether the routine laboratory block improves a coherent probability distribution over survival and three mutually exclusive death groups: UCOD_LEADING=001 diseases of heart; =002 malignant neoplasms; all remaining deaths, including 003-010 and missing cause, as other_or_unknown. Unknown deaths remain events. Heart diseases are NOT synonymous with ischemic heart disease, all cardiovascular disease, or strictly noninfectious death. The residual includes infections and external causes. No WHO population fractions will be used to split individual risk. Detailed ten-category source counts will be retained separately.

## Fixed populations and horizons
Development: NHANES2005-2006 and2007-2008, age40-79, linkage-eligible, positive MEC weights, valid sex and examination follow-up, complete clinical AND routine laboratory predictors under iteration08's fixed definitions. Same participants for clinical and laboratory models. No protein-subset selection and no missing-value imputation. New test: NHANES2009-2010 under identical eligibility; original questionnaire and assay documentation must be checked before interpreting outcomes. No evaluation data are used for feature/penalty selection or calibration.

Primary horizon5years; secondary1 and8years. All surviving test participants must have follow-up through8years for binary evaluation; otherwise stop that horizon and document an amendment, do not silently drop early censored participants. Deaths after a horizon count as event-free for that horizon. Observed zero-month times, if any, use0.5month and are counted. No10/20year extrapolation.

## Models, preprocessing and training
Two panel specifications: the M0c clinical feature list and M1 routine feature list from the frozen model0.6 (SHA256 7b287f3f4778b2127f2226f79d2e503e1e39087306a172e5f524a561964b07d0). This is NEW fitting of cause-specific coefficients, not reuse of published HRs and not a modification of Workbench's default. Clinical comparator was selected in earlier work; selection here precedes the new test.

For each of three causes fit a survey-weighted ridge piecewise exponential model with common slope across intervals0-1,1-5,5-8years and separate unpenalized baseline hazards. All participants contribute at-risk time up to first death/censoring/8years. Training weights normalized to mean1. Missing/cause-unknown deaths are not removed. Optimization uses the full cause-specific likelihood, with ridge(lambda/2)*sum(beta^2); no synthetic deaths or continuity counts. Nonconvergence and zero event counts in required cause/interval cells are explicit blockers, not silent numerical repairs.

Training-only preprocessing: each feature winsorized at unweighted0.5/99.5 percentiles then standardized by its training mean and population SD; zero-variance features retain scale1. Fit preprocessing separately inside validation folds. Penalty candidates1,10,100, selected separately for each panel by mean5year multiclass Brier over two leave-one-survey-cycle-out DEVELOPMENT folds; ties favour the stronger penalty. Then refit selected panel on both development cycles. No after-the-fact marker selection. Coefficients and preprocessing are frozen and hashed before test prediction.

## Coherent probabilities and performance
Within each interval, integrate all cause hazards together: dF_k=S_previous*(1-exp(-sum(lambda)*dt))*lambda_k/sum(lambda). Update survival by exp(-sum(lambda)*dt). Check S+sum(F_k)=1 and monotonicity for every record. Report cause-specific CIF, all-cause risk=sum(F_k), and survival; never independently sum single-cause1-exp(-H_k).

Primary endpoint: weighted multiclass Brier=sum over four outcome classes of squared probability error; report both panels and M1-minus-clinical difference. Secondary: cause-specific and all-cause binary Brier, AUC(cause cases versus everyone else, including competing deaths), observed/expected, calibration bins and mean risks. Do not call AUC diagnostic accuracy or causal evidence. Report all outcomes including unfavourable differences. Sex and age40-59/60-79 breakdowns are descriptive, no subgroup-driven model revision.

Uncertainty:1000 paired rescaled-PSU bootstrap replicates within test survey strata, m-1draws and m/(m-1)rescaling, seed20260911. Percentile intervals conditional on the trained models; training/selection uncertainty is excluded. Zero-event resamples yield undefined AUC recorded and counted, not fabricated. Primary multiclass metric always retained. No recalibration on the test cohort.

## Source and delivery controls
Use official CDC inputs and unit/method documentation. Preserve source checksums and join/cardinality audits; independent mortality reader comparison. Detailed code groups outside target natural-noninfectious scope remain distinct; no claims of solving fine-grained ICD specificity or RU/DE calibration.

Export only code, aggregate metrics, coefficient bundle, audit, and an optional strict research replay; no participant records or per-person predictions in git or public archives. Existing all-cause results/abstract retained as historical; any new abstract revision must describe exact sample, endpoints and estimates. Changed scientific choices after new outcomes require an appended amendment and another validation sample for confirmatory claims.
