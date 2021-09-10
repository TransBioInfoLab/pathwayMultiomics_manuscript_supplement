# AD Multi-Omics with mitch::
# Gabriel Odom
# 2021-08-21

library(tidyverse)
library(mitch)



######  Single-Gene Data  #####################################################

###  SNP  ###
snp_df <- read_csv(
	"../../GWAS_pathway_method/single_gene_SNP_min_pValue_spline_20210823.csv"
) %>%
	select(EnsID = ensembl_gene_id, snp_score = scoreSplineResid)


###  DNAm  ###
dnaM_df <- read_csv(
	"DNAm_meta/single_gene_spline_pValues_20210821.csv"
) %>%
	select(EnsID, dnaM_score = scoreSplineResid)


###  RNAseq  ###
rnaSeq_df <- read_csv(
	"RNASeq_ROSMAP/ROSMAP_only/rnaSeq_single_gene_pVals_20210810.csv"
) %>%
	mutate(rnaSeq_score = -log10(stage_pVal)) %>%
	select(EnsID = ensembl_gene_id, rnaSeq_score)



###  Join  ###
martMatch_df <-
	TCGAbiolinks::get.GRCh.bioMart("hg19") %>%
	as_tibble() %>%
	select(
		EnsID  = ensembl_gene_id,
		symbol = external_gene_name,
		Entrez = entrezgene_id
	) %>%
	distinct()

adMultiOmics_df <-
	snp_df %>%
	full_join(rnaSeq_df, by = "EnsID") %>%
	full_join(dnaM_df, by = "EnsID") %>%
	left_join(martMatch_df, by = "EnsID") %>%
	filter(!is.na(Entrez))

write_csv(
	adMultiOmics_df,
	file = "./AD_single_gene_pValue_scores_Entrez_20210823.csv"
)

rm(list = ls())



######  Use mitch::  ##########################################################
geneSets_ls <- gmt_import("../../gene_sets/c2_cp_entrez_trimmed_20210709.gmt")

adMultiOmics_df <- read_csv(
	file = "./AD_single_gene_pValue_scores_Entrez_20210823.csv"
) %>%
	distinct(Entrez, .keep_all = TRUE) %>%
	column_to_rownames("Entrez") %>%
	select(-EnsID, -symbol)

mitch_out <- mitch_calc(
	x = adMultiOmics_df,
	genesets = geneSets_ls,
	minsetsize = 5L,
	priority = "significance"
)

write_csv(
	mitch_out$enrichment_result,
	file = "./AD_pathway_mitch_results_20210823.csv"
)
