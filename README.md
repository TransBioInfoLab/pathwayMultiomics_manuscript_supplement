# Supplemental Material for pathwayMultiomics Manuscript  
Repository for supplemental material of the pathwayPCA MultiOmics manuscript  
Gabriel Odom and Lily Wang  
2021-08-27  


# Overview of Repository Contents
The first file we include is the final variant of the C2 CP collection we used for the multi-omics analysis: `c2_cp_entrez_trimmed_20210709.gmt`. It includes only gene sets with between 5 and 200 genes.


## Single-Omics: DNA Methylation
The `analysis_DNAm/` folder contains:

- `input_data/`: 
    + `"meta_analysis_no_crossHyb_smoking_ov_comb_p_singleCpG_20210803.csv"`: statistical modelling results for 742 co-methlyated genomic regions
    + `"meta_analysis_single_cpg_sig_no_crossHyb_smoking_with_state_greatAnnot_df.csv"`: statistical modelling results for 3751 individual methylation probes
- `"DNAm_gene_set_enrichment_20210729.R"`: analysis script to perform gene set methylation analysis via the `missMethyl::` package
- `"c2cp_trimmed_missMethyl_pValues_20210803.csv"`: pathway analysis results 
- `"single_gene_spline_pValues_20210821.csv"`: single gene results used for `mitch::` input


## Single-Omics: RNAseq
The `analysis_RNAseq/` folder contains:

- `"rnaSeq_single_gene_pVals_20210810.csv"`: results from single-gene testing
- `"rnaSeq_single_gene_testing_20210810.R"`: script to wrangle the raw data, perform the single gene testing, and then perform the subsequent `fgsea::` analysis
- `"c2cp_trimmed_fgsea_pValues_20210810.csv"`: pathway-level results returned by the `fgsea::` package


## Single-Omics: SNP
The `analysis_SNP/` folder contains: 

- `example_input_data/`: a folder with 10 Excel files. These files are annotated SNP activity within a single pathway. The full set of data files is 2,833 such tables (one for each pathway in the C2 CP collection). We do not include all these files due to repository size restrictions. The full set of 2833 Excel tables would be fed into SAS via the analysis scripts in `src/sas/`.
- `"single_gene_SNP_min_pValue_spline_20210823.csv"`: single gene results used for `mitch::` input
- `src/`: script files. There are 6 SAS scripts to perform the SNP analysis (in `sas/`) and 2 R scripts to wrangle the SAS output:
    + `"MOD_7-10-2021-clean.sas"`: script for pathway analysis of genetic variants
    + `"p-values-estimation_use-empirical-null-clean.R"`: script used to apply the Bacon significance correction
    + `"wrangle_SAS_output_20210720.R"`: script to clean up the results
- `"all_results_bacon_correction_7-21-2021.csv"`: the SNP pathway analysis results after correcting the *p*-values
- `"all_bacon_results_wrangled_20210809.csv"`: SNP pathway analysis results from the SAS macro


## Multi-Omics
The `multiomics/` folder contains the scripts and results for the MiniMax statistic (`MiniMax/`) and for the `mitch::` package (`mitch/`).


