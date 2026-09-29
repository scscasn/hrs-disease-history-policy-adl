# Variable Dictionary

The code reads only the RAND HRS variables required for the stated design. Identifiers are used transiently only to apply the documented duplicate-record alert and are not retained in any public artifact.

| Construct | Source variables used by code | Release role |
|---|---|---|
| Study key | `HHIDPN` | Transient duplicate-alert check only; never released |
| Eligibility and follow-up | `IWSTAT`, age, baseline weight, survey stratum and half-sample | Cohort definition and survey design |
| Baseline function | Five ADL item recodes | Baseline exclusion |
| Follow-up function | Five ADL item recodes | Wave-specific outcome classification |
| Disease histories | Eight `Ever`, flag, raw and question fields | Final and cutoff-policy condition counts |
| Covariates | Sex, race/ethnicity, education, marital status, proxy interview, body mass index | Observation and outcome-model adjustment |

Field prefixes and wave numbers are specified directly in `analysis/build_hrs_analysis_data.py`. This document intentionally omits all data values and record-level outputs.
