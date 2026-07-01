# Step 3 classifier selection and gene-coverage decision v1

Run timestamp: 2026-06-29 22:09:45 IST

## Objective

This step converts the project from ordinary signature construction to an external-validation and transferability benchmark. The goal is to identify published breast-cancer transcriptomic classifiers or mechanistic scores that are sufficiently transparent to audit across the Step 2 expression matrices.

## Main decision

Use a two-tier classifier set:

1. Main computable benchmark set: PAM50, MammaPrint/70-gene, OncotypeDX 21-gene set, Genomic Grade Index, Wang 76-gene signature, EndoPredict, IFN-gamma 18-gene immune profile, and cytolytic activity score.
2. Literature-anchor but not-yet-computable set: Hess pCR 30-gene predictor, Hatzis pCR/survival predictor, HER2DX, and TNBCtype. These are important for the Introduction/Discussion and for showing that pCR prediction has been attempted, but they should not be scored in the main analysis until their exact public formula, feature preprocessing, and cutoffs are recovered.

## Why this is not a new-signature project

No new pCR signature is selected or optimized here. Step 3 only asks whether published classifiers can be transferred across the already-audited public NAC cohorts. This protects the paper frame: external validation, calibration, transferability, and failure modes.

## Immediate Step 4 consequence

Before any AUC/calibration analysis, lock an outcome-blind preprocessing rule per platform: probe-to-gene collapsing, log2 handling, z-score/quantile standardization choice, and missing-gene rule. Do not tune these choices against pCR/RD in validation cohorts.

## Coverage summary

Key: <classifier_id, classifier_name>
          classifier_id                               classifier_name
                 <char>                                        <char>
1: CYTOLYTIC_ACTIVITY_2                      Cytolytic activity score
2:       ENDOPREDICT_11                       EndoPredict 11-gene set
3:               GENE76     Wang 76-gene distant-metastasis signature
4:        GGI_128_PROBE  Genomic Grade Index implementation in genefu
5:              IFNG_18      IFN-gamma-related 18-gene immune profile
6:        MAMMAPRINT_70                MammaPrint / 70-gene signature
7:       ONCOTYPE_DX_21 Oncotype DX 21-gene recurrence score gene set
8:                PAM50            PAM50 intrinsic subtype classifier
   n_matrices_audited min_coverage_pct median_coverage_pct n_matrices_excellent
                <int>            <num>               <num>                <int>
1:                 12                0               100.0                   10
2:                 12                0               100.0                    9
3:                 12                0                84.5                    0
4:                 12                0                90.9                    0
5:                 12                0               100.0                    7
6:                 12                0                82.4                    0
7:                 12                0               100.0                   11
8:                 12                0                96.0                    7
   n_matrices_usable_or_better
                         <int>
1:                          10
2:                          10
3:                          10
4:                          10
5:                          10
6:                           7
7:                          11
8:                          11
                                                                  matrices_below_80pct
                                                                                <char>
1:                                               GSE106977_GPL17586;GSE109710_GPL24546
2:                                               GSE106977_GPL17586;GSE109710_GPL24546
3:                                               GSE106977_GPL17586;GSE109710_GPL24546
4:                                               GSE106977_GPL17586;GSE109710_GPL24546
5:                                               GSE106977_GPL17586;GSE109710_GPL24546
6: GSE106977_GPL17586;GSE109710_GPL24546;GSE20271_GPL96;GSE25066_GPL96;GSE41998_GPL571
7:                                                                  GSE106977_GPL17586
8:                                                                  GSE106977_GPL17586
   n_source_rows_total n_source_rows_with_harmonized_gene_symbol
                 <int>                                     <int>
1:                   2                                         2
2:                  11                                        11
3:                  76                                        61
4:                 128                                       112
5:                  18                                        18
6:                  70                                        55
7:                  21                                        21
8:                  50                                        50
   n_unique_harmonized_genes
                       <int>
1:                         2
2:                        11
3:                        58
4:                        99
5:                        18
6:                        51
7:                        21
8:                        50

## Key reproducibility caveats

- Commercial or clinical-use assays may have public gene lists but non-public weights, calibration, or cutoff rules. Treat those as gene-set or approximate-score benchmarks unless the exact algorithm is public.
- Legacy signatures contain unmapped probes/contigs or old gene symbols. Coverage percentages are calculated against harmonized gene symbols available in Step 2, not against original probe-level assay content.
- Hess and Hatzis are directly relevant pCR papers, but Step 3 v1 does not reconstruct them because the exact scoring rule was not recovered into an auditable local gene list during this pass.
- GSE106977 and GSE109710 are expected stress points because Step 2 already found limited gene space / panel constraints.

## Output files

- classifier_dictionary_v1.csv
- classifier_gene_lists_v1.csv
- classifier_gene_coverage_v1.csv
- classifier_coverage_summary_v1.csv
- classifier_coverage_by_analysis_set_v1.csv
- classifier_reproducibility_risk_v1.csv
- build_classifier_coverage_audit_v1.R
- R_raw_run_log_step3_v1.txt
- R_session_info_step3_v1.txt
