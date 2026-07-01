# breast-nac-transcriptomic-classifier-benchmark
Code and supplementary materials for benchmarking published transcriptomic classifiers for breast cancer NAC pCR prediction.
# External transferability benchmark of transcriptomic classifiers for breast NAC pCR

This repository contains the code and lightweight reproducibility documentation for the manuscript:

**External transferability and failure-mode benchmarking of published transcriptomic classifiers for neoadjuvant breast cancer pathological complete response**

The analysis benchmarks published transcriptomic classifiers or scores across public pretreatment neoadjuvant chemotherapy (NAC) breast cancer cohorts for pathological complete response (pCR) versus residual disease (RD). The repository is organized to match the code used to generate the manuscript tables, figures, and claim-evidence summaries. It does not introduce a new model-building workflow.

## Scope of this release

Included:

- Scripted cohort manifest construction, expression audit, classifier reconstruction audit, classifier scoring, external validation, sensitivity summaries, and manuscript table generation.
- Lightweight CSV, TXT, and Markdown outputs needed to inspect the analysis record.
- Run order, file manifest, data availability, code availability, and reproducibility notes for Journal of Translational Medicine submission.

Not included:

- GEO raw/download-cache expression matrices. These are public and can be retrieved from GEO using the documented accessions and script URLs.
- Word manuscripts, PowerPoint assembly files, and submission artwork source files.

Included as a reproducibility input:

- `manifests/step 3/genefu_2.44.0.tar.gz`, the public package source tarball used by the Step 3 classifier reconstruction audit.

## License and data reuse

The analysis code and repository documentation are released under the MIT License; see `LICENSE`.

The public source datasets remain attributable to NCBI GEO and the original cohort publications. Reuse of the public datasets should cite the relevant GEO accessions and source publications listed in Supplementary Table S1. This repository license does not transfer copyright in the manuscript text, journal article, or original third-party datasets.

## Public data sources

All cohorts are public GEO breast cancer NAC transcriptomic resources. The analysis uses:

GSE25066, GSE20271, GSE32646, GSE41998, GSE50948, GSE66305, GSE163882, GSE194040, GSE106977, GSE109710, and GSE130786.

See `DATA_AVAILABILITY_JTM.md`, `manifests/table1/table1/Table_1_source_audit_v1.csv`, and Supplementary Table S1 under `manifests/final_article_tables` for dataset-level provenance and GEO-linked source publication/dataset citations.

## Analysis boundary

This project does not train a new pCR signature. It reconstructs or approximates published transcriptomic classifiers/scores where public gene lists and scoring rules allow reproducible implementation, then evaluates transferability, calibration, sensitivity behavior, and failure modes across predefined public cohorts.

Key leakage-control rules:

- GSE25066 is used as the discovery/model-reconstruction base.
- External validation cohorts are not used for feature selection, score refitting, cutoff optimization, direction flipping, calibration fitting, or model tuning.
- GSE194040 is retained as a large platform-specific stress test.
- Subtype-only or limited-panel datasets are treated as exploratory context, not primary validation evidence.

## Software environment

The manuscript-generating run logs record:

- R 4.5.0 on Windows 10 x64.
- Main R packages: `data.table` 1.17.0 and `ggplot2` 3.5.2.
- Optional table export packages: `openxlsx` or `writexl`, if available.
- Step 1 manifest construction uses Windows PowerShell.

See `docs/R_REQUIREMENTS.md` for details.

## Quick run order

Run commands from the repository root after public data access is available:

```powershell
powershell -ExecutionPolicy Bypass -File code/scripts/build_pcr_manifest.ps1 -OutDir "manifests/step 1"
Rscript code/scripts/build_expression_audit_v1.R "manifests/step 1" "manifests/step 2"
Rscript code/scripts/build_classifier_coverage_audit_v1.R .
Rscript code/scripts/build_classifier_scores_step4_v1.R .
Rscript code/scripts/build_external_validation_step5_v1.R .
Rscript code/scripts/build_sensitivity_step6_v1.R .
Rscript code/manuscript_tables/build_table1_public_nac_cohorts_v2.R
Rscript code/manuscript_tables/build_final_article_tables_v1.R
```

For explanations of each step and expected outputs, see `RUN_ORDER.md`.

## Main locked results

The locked primary ranking is stored in:

- `manifests/step 6/primary_result_lock_v1.csv`
- `manifests/step 6/final_classifier_ranking_locked_v1.csv`

The top three locked primary median external AUC values are:

1. IFN-gamma 18-gene mean z score: median external AUC 0.678.
2. PAM50 non-luminal affinity score: median external AUC 0.670.
3. OncotypeDX approximate weighted expression score: median external AUC 0.657.

These results should be interpreted as modest external transferability signals, not as evidence of clinical utility.

## Journal-specific note

Journal of Translational Medicine asks authors to provide datasets supporting conclusions in repositories or supporting files, include data availability information, fully reference public datasets, and provide software details including project name, project home page, archived version/identifier, operating system, programming language, requirements, license, and use restrictions. See:

https://link.springer.com/journal/12967/submission-guidelines

This package includes draft text for the manuscript sections in:

- `DATA_AVAILABILITY_JTM.md`
- `CODE_AVAILABILITY_JTM.md`

Before submission, replace placeholder fields for the GitHub repository URL and archived Zenodo DOI. The code license is MIT License.
