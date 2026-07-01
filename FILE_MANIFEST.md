# File manifest

This manifest maps the manuscript-generating files included in this GitHub release package.

## Code files

| Release path | Original project file | Purpose |
|---|---|---|
| `code/scripts/build_pcr_manifest.ps1` | `scripts/build_pcr_manifest.ps1` | Build public pCR/RD phenotype manifest from GEO metadata |
| `code/scripts/build_expression_audit_v1.R` | `manifests/step 2/build_expression_audit_v1.R` | Audit public expression matrices and feature annotation |
| `code/scripts/build_classifier_coverage_audit_v1.R` | `manifests/step 3/build_classifier_coverage_audit_v1.R` | Audit classifier gene lists, reconstructability, and gene coverage |
| `code/scripts/inspect_genefu_signatures_v1.R` | `manifests/step 3/inspect_genefu_signatures_v1.R` | Inspect genefu signature objects |
| `code/scripts/build_classifier_scores_step4_v1.R` | `manifests/step 4/build_classifier_scores_step4_v1.R` | Generate outcome-blind classifier scores |
| `code/scripts/build_external_validation_step5_v1.R` | `manifests/step 5/build_external_validation_step5_v1.R` | External validation, calibration, failure-mode tables, and Step 5 figures |
| `code/scripts/build_sensitivity_step6_v1.R` | `manifests/step 6/build_sensitivity_step6_v1.R` | Sensitivity summaries, locked ranking, calibration failure summary, and claim-evidence table |
| `code/manuscript_tables/build_table1_public_nac_cohorts_v2.R` | `manifests/table1/build_table1_public_nac_cohorts_v2.R` | Build manuscript Table 1 |
| `code/manuscript_tables/build_final_article_tables_v1.R` | `manifests/final_article_tables/build_final_article_tables_v1.R` | Build manuscript Tables 2-4 and supplementary tables |
| `code/manuscript_tables/qc_user_arranged_figures_tables_v1.R` | `manifests/final_article_tables/qc_user_arranged_figures_tables_v1.R` | QC table and figure package consistency |

## Lightweight reproducibility outputs

The release includes CSV/TXT/MD outputs from the original manuscript-generating directories:

- `manifests/step 1`
- `manifests/step 2`
- `manifests/step 3`
- `manifests/step 4`
- `manifests/step 5`
- `manifests/step 6`
- `manifests/table1`
- `manifests/final_article_tables`

The original directory names, including spaces in `step 1`, `step 2`, and later step folders, are intentionally retained because the scripts read and write those paths.

Supplementary Table S1 in `manifests/final_article_tables` includes the final `Source publication / dataset citation` column with GEO-linked PMIDs/DOIs for the public datasets.

Excluded file classes:

- GEO download-cache files under `manifests/step 4/download_cache`.
- Duplicate desktop table-package copies of the final article tables.
- Word manuscripts and OOXML build folders.
- PowerPoint files and editable manuscript figure decks.
- Submission artwork source files; code-generated Step 5/6 figures may be retained as analysis-output examples.
- Failed-attempt logs unless specifically needed for provenance.

## Archived release DOI or persistent identifier: to be added after archival.

- GitHub repository URL.
- Archived release DOI or persistent identifier.

The software license has been selected: MIT License, provided in `LICENSE`.
