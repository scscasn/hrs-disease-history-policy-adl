# HRS Disease History Policy and ADL Association

This repository contains code and documentation for the manuscript **Retrospective Chronic-Condition Revisions and Multimorbidity–ADL Associations in the Health and Retirement Study: A Cohort Study of Practical Equivalence Among Older Survivors**.

## Scope

The package preserves the release-version analytic specification, code, correction record, and reproducibility instructions. It does not include Health and Retirement Study (HRS) microdata, row-level derived datasets, identifiers, credentials, or restricted materials.

## Data access and permitted use

Reproduction requires an independently obtained, authorized copy of the RAND HRS Longitudinal File 1992–2022 v1 Early Release data product, the applicable official HRS data-correction alert, and adherence to all HRS Conditions of Use. Users must provide local source and correction-list paths to `analysis/build_hrs_analysis_data.py` themselves. HRS data, row-level derivatives, and record-specific identifiers must not be uploaded, shared, or committed to this repository.

The release-version analysis was manually rerun by the authors in a local non-AI environment. This package should likewise be run locally without providing HRS person-level data to AI systems or allowing AI tools to index directories containing HRS person-level data.

## Contents

- `analysis/`: protected-local data construction, primary analysis, independent reproduction, and table/figure scripts.
- `documentation/`: frozen analysis specification and reproduction record.
- `metadata/`: software requirements and a disclosure-oriented variable dictionary.

## Reproduction sequence

1. Obtain authorized RAND HRS 1992–2022 v1 Early Release data and verify the documented SHA-256 locally.
2. In a local protected workspace, create an untracked correction-list file from the applicable official HRS data alert, then run `python build_hrs_analysis_data.py --source /path/to/randhrs1992_2022v1_STATA.zip --official-duplicate-id-file /path/to/local_correction_list.txt --output-dir derived`.
3. Run `Rscript run_policy_equivalence_analysis.R /path/to/run_directory`.
4. Run `Rscript audit_primary_independent.R /path/to/run_directory` and compare the reproduced estimate with the declared tolerance.
5. Run `Rscript build_tables_and_figures.R /path/to/run_directory` for manuscript displays.

## Version

Version 1.0.0 corresponds to the frozen analysis and manuscript package dated 29 September 2026. The public release deliberately omits record-specific correction identifiers; authorized users must obtain them from the applicable HRS data alert and keep them local. The code-only release is available at https://github.com/scscasn/hrs-disease-history-policy-adl under the public release tag `v1.0.0`.
