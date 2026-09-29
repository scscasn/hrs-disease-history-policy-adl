# Proposed Q1 SAP lock v2 — HRS information cutoff and subsequent functional status

Version: `HRS-ASOF-ADL-SAP-v2.1`  
Status: proposed, result-blind, not yet user-approved or executable.  
Frozen before: any eligible-sample disease-policy cross-tabulation, disease–ADL contrast, model coefficient, discrimination/calibration comparison or p-value.

## Scientific question and claim boundary

Question: Does a chronic-condition information policy restricted to what is identifiable by the HRS index wave materially change the association between multimorbidity and newly observed five-task ADL difficulty among follow-up survivors, compared with the final release's retrospectively revised `Ever` indicators?

Claim type: methodological epidemiology and measurement sensitivity. No causal, treatment, clinical-diagnostic, exact-onset, permanent-disability, deployable-prediction or historical-download claim is allowed. A coding-policy comparison evaluates association sensitivity; it does not establish which policy is clinically true.

## Population, time zero and analysis pairs

- Source: read-only RAND HRS 1992–2022 v1.
- Primary pair: Wave 14 (2018) baseline → Wave 15 (2020) follow-up.
- Time-horizon sensitivity: Wave 15 (2020) → Wave 16 (2022). It is not called an independent replication because the final release provides different later-revision opportunities and many people can appear in both pairs.
- Baseline age ≥65 years; all five item-level ADL-difficulty recodes observed as 0/1 and all zero; positive combined community/nursing-home respondent weight; one record per `HHIDPN` per pair.
- Structural primary envelope: 7,046 baseline records. The later-pair envelope is 7,141. These counts precede the exposure-ambiguity exclusion.

## Frozen disease-information policies

Eight conditions: hypertension, diabetes, cancer, lung disease, heart disease, stroke, psychiatric problems and arthritis.

1. **Final retrospective policy:** use the final release's same-wave `RwCONDITIONE` values.
2. **Partially identified cutoff policy:** set `F=6` to 1; set `F=4/5` to unknown; otherwise retain observed final `Ever` 0/1.

The primary target is explicitly restricted to baseline-eligible people who survive to follow-up and whose eight cutoff states are identifiable. Because `F=4/5` is learned from later information, this exclusion is itself future-dependent; results cannot be extrapolated to every baseline respondent or described as a complete historical cutoff reconstruction. Both policies use exactly the same people, survey/observation weights and standardization distribution.

Ambiguous `F=4/5` records are reintroduced in two prespecified **extreme-assignment scenarios**, first assigning every ambiguous condition 0 and then 1. These scenarios test sensitivity to the primary exclusion, but they are not claimed to bound a nonlinear regression coefficient or the policy contrast: mixed assignments could be more extreme. No data-driven condition dropping is allowed.

For each policy, form a count from 0 to 8 and fixed categories 0, 1, 2 and ≥3. The continuous count is primary; category contrasts are secondary because the ≥3 group can be sparse after exclusions.

## Outcome and estimand

The target population for the primary functional estimand is baseline-eligible people who survive to the follow-up wave.

- Positive follow-up functional status: at least one observed item-level ADL difficulty equals 1, even if another item is missing.
- Negative status: all five items observed and all equal 0.
- Otherwise: functional status unknown.

For each policy, fit the same modified-Poisson model among observed follow-up survivors:

`log E(Y) = α + βM + γ1(age_c/10) + γ2(age_c/10)^2 + sex + race/ethnicity + education_years + marital/partner status + baseline_proxy`,

where `age_c = age − 75` and `M` is the policy-specific 0–8 condition count, entered linearly. The standardization population is the same exposure-complete survivor target, using its baseline survey-weighted covariate distribution.

The **single primary estimand** is `ΔlogPR = β_cutoff − β_final`, the difference in log prevalence ratios per additional condition. Report `ΔlogPR`, its 95% interval and `exp(ΔlogPR)` (the ratio of the two per-condition PRs). No alternative primary ratio/difference will be selected after results are seen. The linear-per-condition assumption is checked secondarily with the frozen 0/1/2/≥3 categories, not used to replace the primary estimand.

A prespecified practical-equivalence region is `|ΔlogPR| < log(1.10)`, corresponding to less than a 10% relative difference between the two per-condition PRs. A 95% interval wholly inside this region is an informative near-null; a point estimate outside it is potentially consequential but still requires uncertainty and sensitivity review; an interval spanning the region is inconclusive. This threshold is an interpretation aid fixed before analysis, not a publication guarantee. The functional outcome is not called cumulative incidence or a first event.

Key secondary estimands are standardized prevalence differences for categories 1, 2 and ≥3 versus 0; the adverse follow-up composite of death or known ADL difficulty in the full baseline population; and death alone as a separate descriptive competing state.

The result-blind margins are: Wave 14→15, 640 known-positive functional states, 5,378 known-negative, 397 deaths and 631 alive with unknown status/nonresponse; Wave 15→16, 636 known-positive, 5,046 known-negative, 425 deaths and 1,033 alive with unknown status/nonresponse. These are unstratified outcome/observation margins, not policy effects.

## Observation, missingness and covariates

Among survivors, one common set of stabilized inverse-probability-of-observation weights will be estimated using index-available age, sex, race/ethnicity, education, marital/partner status, baseline proxy, BMI and the identified cutoff count. The two policies use these exact same weights, people and target distribution.

- The observation-probability denominator is a logistic model containing `(age−75)/10`, its square, sex, race/ethnicity, education years, marital/partner status, baseline proxy, BMI and the cutoff count. Its numerator is the design-weighted overall probability of known functional status among target survivors.
- Education and BMI are median-imputed within wave with separate missing indicators in the observation model; categorical missingness uses explicit missing levels. The outcome model uses the same fixed education fill/indicator rule. BMI is excluded from the primary outcome model and added only in a secondary model.
- Primary interpretation assumes missing at random conditional on the fixed predictors.
- Stabilized observation weights are truncated at 1st/99th percentiles; no truncation and 2.5th/97.5th are sensitivities.
- A prespecified delta-adjusted pattern-mixture sensitivity includes **all** alive people whose functional state is unknown, including nonresponders and interviewed people whose partial items do not identify positive or negative status. A shared survey-weighted logistic outcome model based on cutoff count and the fixed covariates predicts their outcome probabilities; odds are multiplied by 0.5, 1, 2 and 3. For each multiplier the same completed expected-outcome vector is analyzed under both coding policies, and the multiplier at which the practical conclusion changes is reported. Alive unknown states in the death/difficulty composite are never coded negative. Multiple imputation of baseline covariates is secondary and never presents an imputed exposure or outcome as observed.

## Survey design and inference

Construct the survey design on the full baseline weighted respondent frame using the combined respondent/nursing-home weight, `RAESTRAT` and `RAEHSAMP`, then analyze the age/function cohort as a survey domain. The audited full frames contain all 80 strata and 160 stratum/half-sample units with no singleton strata; singleton units created only by the analytic subdomain must not be treated as if the design had been built after filtering.

Fit the identical survey-weighted modified-Poisson formula under the two policies. Use Taylor-linearized variance for policy-specific estimates. For `ΔlogPR`, convert the full-frame two-unit-per-stratum design to deterministic Fay BRR replicate weights with `rho=0.5` and `mse=TRUE`, then domain-subset. Within every replicate, refit the common observation model, rebuild stabilized observation weights, refit both policy models and restandardize them over the same replicate-weighted target. The two policies always share the same replicate.

The analysis must record software versions, the replicate-weight matrix hash, convergence and every failed fit. No failed replicate is silently replaced. If more than 5% of scheduled replicates fail, or the two policies do not use identical successful replicates, no primary interval is released until code/model repair and repeat review. Statistical code review is mandatory before any effect reaches the publication gate.

## Prespecified sensitivities

1. `F=4/5` all-0 and all-1 extreme-assignment scenarios; no formal-bound claim.
2. Wave 15→16 time-horizon sensitivity, explicitly not an equivalent replication.
3. Self interviews only; all respondents with proxy adjustment.
4. Complete five-item follow-up outcome versus the partial-item positive rule.
5. Unweighted analysis as a design sensitivity.
6. BMI-adjusted secondary outcome model.
7. Leave-one-condition-out and condition-level changes as exploratory analyses with false-discovery control.
8. Interview-mode sensitivity where a documented mode proxy is available; no undocumented cross-wave mode variable is invented.

No AUC, calibration or deployable-prediction claim is planned: the final retrospective policy contains information unavailable at the index and is therefore not a valid deployment comparator.

## Research-value and stop rules

The exact prior search found close HRS adjudication and disability precedents but not this time-cutoff consequence analysis. The paper is viable only if it changes or precisely resolves a real data-construction decision.

- Continue at a JIF≥5 aspirational tier only if the paired consequence crosses the prespecified 10% relative-change threshold with persuasive precision, or the full interval lies inside the equivalence region while extreme-assignment, MNAR and time-horizon sensitivities agree.
- A precise near-null result is publishable evidence if it answers the method-choice question; it is not an automatic failure.
- Downgrade to a specialist ageing/measurement journal when the construct is valid but the consequence is modest.
- Recommend abandonment when intervals cannot distinguish meaningful change, observation-process assumptions dominate, or the conclusion reverses across extreme-assignment scenarios/time horizons.
- Do not select whichever association, subgroup or metric changes most after seeing results.

## Conditional journal ladder

- Aspirational if consequence is strong and robust: *International Journal of Epidemiology*; *Age and Ageing* only when the functional implication is central.
- Realistic specialist fallback for a valid but modest methods lesson: *The Journals of Gerontology: Series B* or a comparable ageing/measurement journal.

This is a proposed Q1 lock, not an acceptance forecast. The true submission interval will be disclosed after Q2 results and formal publication-gate review. No deletion follows automatically: the user first receives the real results and realistic tier, then decides whether to abandon and clean up.
