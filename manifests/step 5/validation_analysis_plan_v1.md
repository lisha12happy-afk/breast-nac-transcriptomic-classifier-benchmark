# Step 5 external validation analysis plan v1

Generated: 2026-06-30 19:29:20 IST

## Objective

Evaluate whether pre-specified transcriptomic classifier scores from Step 4 transfer across independent breast cancer neoadjuvant chemotherapy cohorts for pCR/RD prediction.

## Primary endpoint

`endpoint_binary_pcr1_rd0`, where pCR = 1 and residual disease = 0.

## Primary score rule

Discrimination is reported as discovery-aligned AUC. Direction is learned only in the discovery/calibration base GSE25066; external cohorts are not used to flip score direction.

## Primary analysis set

External core validation cohorts with `analysis_set == core_full_benchmark` and `score_status == scored_preferred_coverage`.

## Secondary/exploratory analysis sets

- GSE194040 platforms are reported as large stress tests.
- Limited gene-space/panel datasets are reported separately.
- Low-coverage scores are flagged as exploratory.

## Metrics

- ROC-AUC with bootstrap 95% CI.
- Average precision / PR-AUC with pCR prevalence as baseline.
- Discovery-trained logistic calibration: Brier score, calibration intercept, calibration slope, observed/expected ratio.
- Failure-mode summaries by platform, analysis set, and score coverage status.

## Leakage control

No external validation outcome is used for direction flipping, score scaling, cutoff selection, calibration fitting, or model tuning.
