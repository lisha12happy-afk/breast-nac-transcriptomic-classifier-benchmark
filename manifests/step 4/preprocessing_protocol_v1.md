# Step 4 outcome-blind preprocessing and classifier scoring protocol v1

Generated: 2026-06-29 23:18:17 IST

## Fixed analysis boundary

- This step does not build a new pCR signature and does not select genes using pCR/RD.
- The purpose is to convert each public expression matrix into deterministic published-classifier or mechanistic scores.
- GSE194040 platforms are kept separate.
- GSE106977 and GSE109710 are retained as limited gene-space/panel stress checks but should not drive main conclusions.

## Expression preprocessing

1. Expression matrices are read from the public URLs already audited in Step 2.
2. Manifest sample matching uses GSM/sample title/patient aliases from Step 1 and the same normalization rule as Step 2.
3. Features are mapped to Step 2 feature-index gene symbols.
4. Probes mapping to multiple genes are excluded from scoring to avoid ambiguous feature attribution.
5. Duplicate single-gene probes are collapsed by the arithmetic mean across probes for each sample. This is deterministic and outcome-blind.
6. Matrices labelled as likely unlogged intensity/count/TPM, or with sampled maxima >50, are transformed as log2(x+1) after truncating negative values to zero. Matrices already on a log-like continuous scale are kept as-is.
7. Gene expression is standardized within each matrix across matched manifest samples: z = (expression - gene mean) / gene SD.
8. Genes with zero or non-finite SD are dropped before scoring.

## Leakage control

- pCR/RD, receptor status, treatment arm, and planned-role labels are never used to choose preprocessing parameters.
- Coverage thresholds are pre-specified in classifier_scoring_spec_v1.csv.
- Low-coverage scores are flagged as exploratory or not scored; they are not silently mixed into the main benchmark.

## Scoring outputs

- `classifier_score_matrix_v1.csv`: long sample-by-classifier score table.
- `pam50_subtype_assignments_v1.csv`: PAM50 centroid correlations and subtype calls.
- `classifier_score_qc_v1.csv`: per-dataset per-classifier coverage and scoring status.
- `dataset_preprocessing_decisions_v1.csv`: scale transform, matched samples, gene-level matrix dimensions.
- `sample_score_match_qc_v1.csv`: manifest-to-expression sample matching status.
