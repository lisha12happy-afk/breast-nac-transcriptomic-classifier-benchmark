# Step 6 sensitivity analysis and primary result lock v1

Generated: 2026-06-30 19:31:14 IST

## Purpose

Lock the main external-validation conclusion and test whether it changes under prespecified sensitivity scenarios.

## Locked primary analysis

- `performance_role == external_core_validation`
- `score_status == scored_preferred_coverage`
- Discovery base GSE25066 is excluded from external validation ranking.
- GSE194040 is reported as stress-test only.
- Low-coverage exploratory scores and limited-panel datasets are not mixed into the primary result.

## Sensitivity scenarios

- Primary external core, preferred coverage: performance_role == external_core_validation; score_status == scored_preferred_coverage
- External core, include low-coverage exploratory: performance_role == external_core_validation; score_status in preferred or low-coverage exploratory
- External core microarray only, preferred coverage: performance_role == external_core_validation; score_status == preferred; platform in GPL96/GPL570/GPL571
- External core excluding GPL96, preferred coverage: performance_role == external_core_validation; score_status == preferred; platform != GPL96
- GSE194040 large stress-test only, preferred coverage: performance_role == large_stress_test; score_status == preferred
- All non-discovery matrices, preferred coverage: performance_role != discovery_calibration_base; score_status == preferred
- All non-discovery matrices, include low coverage: performance_role != discovery_calibration_base; score_status in preferred or low-coverage exploratory

## Interpretation boundary

Step 6 can support a paper claim about transferability, instability, and failure modes. It cannot support a claim of clinical utility, prospective predictive validity, or superiority over commercial assays.
