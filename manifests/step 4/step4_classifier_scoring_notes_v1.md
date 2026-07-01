# Step 4 classifier scoring execution notes v1

Generated: 2026-06-29 23:25:59 IST

## Completed outputs

- preprocessing_protocol_v1.md
- classifier_scoring_spec_v1.csv
- dataset_preprocessing_decisions_v1.csv
- sample_score_match_qc_v1.csv
- gene_level_matrix_qc_v1.csv
- classifier_score_matrix_v1.csv
- classifier_score_qc_v1.csv
- pam50_subtype_assignments_v1.csv
- classifier_score_errors_v1.csv
- build_classifier_scores_step4_v1.R
- R_raw_run_log_step4_v1.txt
- R_session_info_step4_v1.txt

## Important interpretation

These scores are pre-specified transferability benchmarks, not newly trained pCR predictors. Commercial or clinical classifiers are represented only by transparent gene-level approximate scores when the exact clinical algorithm is not public.

## Error summary

No dataset-level scoring errors recorded.

## Score QC summary

    dataset_accession       matrix_label        classifier_id
               <char>             <char>               <char>
 1:          GSE25066     GSE25066_GPL96                PAM50
 2:          GSE25066     GSE25066_GPL96        MAMMAPRINT_70
 3:          GSE25066     GSE25066_GPL96       ONCOTYPE_DX_21
 4:          GSE25066     GSE25066_GPL96        GGI_128_PROBE
 5:          GSE25066     GSE25066_GPL96               GENE76
 6:          GSE25066     GSE25066_GPL96       ENDOPREDICT_11
 7:          GSE25066     GSE25066_GPL96              IFNG_18
 8:          GSE25066     GSE25066_GPL96 CYTOLYTIC_ACTIVITY_2
 9:          GSE20271     GSE20271_GPL96                PAM50
10:          GSE20271     GSE20271_GPL96        MAMMAPRINT_70
11:          GSE20271     GSE20271_GPL96       ONCOTYPE_DX_21
12:          GSE20271     GSE20271_GPL96        GGI_128_PROBE
13:          GSE20271     GSE20271_GPL96               GENE76
14:          GSE20271     GSE20271_GPL96       ENDOPREDICT_11
15:          GSE20271     GSE20271_GPL96              IFNG_18
16:          GSE20271     GSE20271_GPL96 CYTOLYTIC_ACTIVITY_2
17:          GSE32646    GSE32646_GPL570                PAM50
18:          GSE32646    GSE32646_GPL570        MAMMAPRINT_70
19:          GSE32646    GSE32646_GPL570       ONCOTYPE_DX_21
20:          GSE32646    GSE32646_GPL570        GGI_128_PROBE
21:          GSE32646    GSE32646_GPL570               GENE76
22:          GSE32646    GSE32646_GPL570       ENDOPREDICT_11
23:          GSE32646    GSE32646_GPL570              IFNG_18
24:          GSE32646    GSE32646_GPL570 CYTOLYTIC_ACTIVITY_2
25:          GSE41998    GSE41998_GPL571                PAM50
26:          GSE41998    GSE41998_GPL571        MAMMAPRINT_70
27:          GSE41998    GSE41998_GPL571       ONCOTYPE_DX_21
28:          GSE41998    GSE41998_GPL571        GGI_128_PROBE
29:          GSE41998    GSE41998_GPL571               GENE76
30:          GSE41998    GSE41998_GPL571       ENDOPREDICT_11
31:          GSE41998    GSE41998_GPL571              IFNG_18
32:          GSE41998    GSE41998_GPL571 CYTOLYTIC_ACTIVITY_2
33:          GSE50948    GSE50948_GPL570                PAM50
34:          GSE50948    GSE50948_GPL570        MAMMAPRINT_70
35:          GSE50948    GSE50948_GPL570       ONCOTYPE_DX_21
36:          GSE50948    GSE50948_GPL570        GGI_128_PROBE
37:          GSE50948    GSE50948_GPL570               GENE76
38:          GSE50948    GSE50948_GPL570       ENDOPREDICT_11
39:          GSE50948    GSE50948_GPL570              IFNG_18
40:          GSE50948    GSE50948_GPL570 CYTOLYTIC_ACTIVITY_2
41:          GSE66305    GSE66305_GPL570                PAM50
42:          GSE66305    GSE66305_GPL570        MAMMAPRINT_70
43:          GSE66305    GSE66305_GPL570       ONCOTYPE_DX_21
44:          GSE66305    GSE66305_GPL570        GGI_128_PROBE
45:          GSE66305    GSE66305_GPL570               GENE76
46:          GSE66305    GSE66305_GPL570       ENDOPREDICT_11
47:          GSE66305    GSE66305_GPL570              IFNG_18
48:          GSE66305    GSE66305_GPL570 CYTOLYTIC_ACTIVITY_2
49:         GSE163882 GSE163882_GPL18573                PAM50
50:         GSE163882 GSE163882_GPL18573        MAMMAPRINT_70
51:         GSE163882 GSE163882_GPL18573       ONCOTYPE_DX_21
52:         GSE163882 GSE163882_GPL18573        GGI_128_PROBE
53:         GSE163882 GSE163882_GPL18573               GENE76
54:         GSE163882 GSE163882_GPL18573       ENDOPREDICT_11
55:         GSE163882 GSE163882_GPL18573              IFNG_18
56:         GSE163882 GSE163882_GPL18573 CYTOLYTIC_ACTIVITY_2
57:         GSE194040 GSE194040_GPL20078                PAM50
58:         GSE194040 GSE194040_GPL20078        MAMMAPRINT_70
59:         GSE194040 GSE194040_GPL20078       ONCOTYPE_DX_21
60:         GSE194040 GSE194040_GPL20078        GGI_128_PROBE
61:         GSE194040 GSE194040_GPL20078               GENE76
62:         GSE194040 GSE194040_GPL20078       ENDOPREDICT_11
63:         GSE194040 GSE194040_GPL20078              IFNG_18
64:         GSE194040 GSE194040_GPL20078 CYTOLYTIC_ACTIVITY_2
65:         GSE194040 GSE194040_GPL30493                PAM50
66:         GSE194040 GSE194040_GPL30493        MAMMAPRINT_70
67:         GSE194040 GSE194040_GPL30493       ONCOTYPE_DX_21
68:         GSE194040 GSE194040_GPL30493        GGI_128_PROBE
69:         GSE194040 GSE194040_GPL30493               GENE76
70:         GSE194040 GSE194040_GPL30493       ENDOPREDICT_11
71:         GSE194040 GSE194040_GPL30493              IFNG_18
72:         GSE194040 GSE194040_GPL30493 CYTOLYTIC_ACTIVITY_2
73:         GSE106977 GSE106977_GPL17586                PAM50
74:         GSE106977 GSE106977_GPL17586        MAMMAPRINT_70
75:         GSE106977 GSE106977_GPL17586       ONCOTYPE_DX_21
76:         GSE106977 GSE106977_GPL17586        GGI_128_PROBE
77:         GSE106977 GSE106977_GPL17586               GENE76
78:         GSE106977 GSE106977_GPL17586       ENDOPREDICT_11
79:         GSE106977 GSE106977_GPL17586              IFNG_18
80:         GSE106977 GSE106977_GPL17586 CYTOLYTIC_ACTIVITY_2
81:         GSE109710 GSE109710_GPL24546                PAM50
82:         GSE109710 GSE109710_GPL24546        MAMMAPRINT_70
83:         GSE109710 GSE109710_GPL24546       ONCOTYPE_DX_21
84:         GSE109710 GSE109710_GPL24546        GGI_128_PROBE
85:         GSE109710 GSE109710_GPL24546               GENE76
86:         GSE109710 GSE109710_GPL24546       ENDOPREDICT_11
87:         GSE109710 GSE109710_GPL24546              IFNG_18
88:         GSE109710 GSE109710_GPL24546 CYTOLYTIC_ACTIVITY_2
89:         GSE130786  GSE130786_GPL6480                PAM50
90:         GSE130786  GSE130786_GPL6480        MAMMAPRINT_70
91:         GSE130786  GSE130786_GPL6480       ONCOTYPE_DX_21
92:         GSE130786  GSE130786_GPL6480        GGI_128_PROBE
93:         GSE130786  GSE130786_GPL6480               GENE76
94:         GSE130786  GSE130786_GPL6480       ENDOPREDICT_11
95:         GSE130786  GSE130786_GPL6480              IFNG_18
96:         GSE130786  GSE130786_GPL6480 CYTOLYTIC_ACTIVITY_2
    dataset_accession       matrix_label        classifier_id
    n_scores_nonmissing coverage_pct              score_status_values
                  <int>        <num>                           <char>
 1:                 508         84.0        scored_preferred_coverage
 2:                 508         70.6  scored_low_coverage_exploratory
 3:                 508         93.8        scored_preferred_coverage
 4:                 508         87.9        scored_preferred_coverage
 5:                 508         82.8        scored_preferred_coverage
 6:                 508        100.0        scored_preferred_coverage
 7:                 508         88.9        scored_preferred_coverage
 8:                 508        100.0        scored_preferred_coverage
 9:                 178         84.0        scored_preferred_coverage
10:                 178         70.6  scored_low_coverage_exploratory
11:                 178         93.8        scored_preferred_coverage
12:                 178         87.9        scored_preferred_coverage
13:                 178         82.8        scored_preferred_coverage
14:                 178        100.0        scored_preferred_coverage
15:                 178         88.9        scored_preferred_coverage
16:                 178        100.0        scored_preferred_coverage
17:                 115         94.0        scored_preferred_coverage
18:                 115         80.4        scored_preferred_coverage
19:                 115         93.8        scored_preferred_coverage
20:                 115         87.9        scored_preferred_coverage
21:                 115         82.8        scored_preferred_coverage
22:                 115        100.0        scored_preferred_coverage
23:                 115        100.0        scored_preferred_coverage
24:                 115        100.0        scored_preferred_coverage
25:                 279         84.0        scored_preferred_coverage
26:                 279         70.6  scored_low_coverage_exploratory
27:                 279         93.8        scored_preferred_coverage
28:                 279         87.9        scored_preferred_coverage
29:                 279         82.8        scored_preferred_coverage
30:                 279        100.0        scored_preferred_coverage
31:                 279         88.9        scored_preferred_coverage
32:                 279        100.0        scored_preferred_coverage
33:                 156         94.0        scored_preferred_coverage
34:                 156         80.4        scored_preferred_coverage
35:                 156         93.8        scored_preferred_coverage
36:                 156         87.9        scored_preferred_coverage
37:                 156         82.8        scored_preferred_coverage
38:                 156        100.0        scored_preferred_coverage
39:                 156        100.0        scored_preferred_coverage
40:                 156        100.0        scored_preferred_coverage
41:                  88         94.0        scored_preferred_coverage
42:                  88         80.4        scored_preferred_coverage
43:                  88         93.8        scored_preferred_coverage
44:                  88         87.9        scored_preferred_coverage
45:                  88         82.8        scored_preferred_coverage
46:                  88        100.0        scored_preferred_coverage
47:                  88        100.0        scored_preferred_coverage
48:                  88        100.0        scored_preferred_coverage
49:                 222         96.0        scored_preferred_coverage
50:                 222         84.3        scored_preferred_coverage
51:                 222        100.0        scored_preferred_coverage
52:                 222         90.9        scored_preferred_coverage
53:                 222         84.5        scored_preferred_coverage
54:                 222        100.0        scored_preferred_coverage
55:                 222        100.0        scored_preferred_coverage
56:                 222        100.0        scored_preferred_coverage
57:                 654         96.0        scored_preferred_coverage
58:                 651         84.3        scored_preferred_coverage
59:                 654        100.0        scored_preferred_coverage
60:                 654         90.9        scored_preferred_coverage
61:                 654         84.5        scored_preferred_coverage
62:                 654        100.0        scored_preferred_coverage
63:                 654        100.0        scored_preferred_coverage
64:                 654        100.0        scored_preferred_coverage
65:                 334         96.0        scored_preferred_coverage
66:                 328         84.3        scored_preferred_coverage
67:                 334        100.0        scored_preferred_coverage
68:                 334         90.9        scored_preferred_coverage
69:                 282         82.8        scored_preferred_coverage
70:                 334        100.0        scored_preferred_coverage
71:                 334        100.0        scored_preferred_coverage
72:                 334        100.0        scored_preferred_coverage
73:                   0          0.0 not_scored_insufficient_coverage
74:                   0          0.0 not_scored_insufficient_coverage
75:                   0          0.0 not_scored_insufficient_coverage
76:                   0          0.0 not_scored_insufficient_coverage
77:                   0          0.0 not_scored_insufficient_coverage
78:                   0          0.0 not_scored_insufficient_coverage
79:                   0          0.0 not_scored_insufficient_coverage
80:                   0          0.0 not_scored_insufficient_coverage
81:                 173        100.0        scored_preferred_coverage
82:                   0         21.6 not_scored_insufficient_coverage
83:                 173        100.0        scored_preferred_coverage
84:                   0         27.3 not_scored_insufficient_coverage
85:                   0         12.1 not_scored_insufficient_coverage
86:                   0         25.0 not_scored_insufficient_coverage
87:                   0         44.4 not_scored_insufficient_coverage
88:                   0          0.0 not_scored_insufficient_coverage
89:                 110         94.0        scored_preferred_coverage
90:                 107         82.4        scored_preferred_coverage
91:                 108        100.0        scored_preferred_coverage
92:                 110         89.9        scored_preferred_coverage
93:                 105         82.8        scored_preferred_coverage
94:                 107         87.5        scored_preferred_coverage
95:                 110        100.0        scored_preferred_coverage
96:                 110        100.0        scored_preferred_coverage
    n_scores_nonmissing coverage_pct              score_status_values
