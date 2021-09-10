# AD Multi-Omics
# Gabriel Odom
# 2021-08-09

# We have extracted pathway-level p-values for SNP (GWAS), DNA methylation, and
#   RNAseq for MSigDB's C2 CP collection



######  Data  #################################################################
library(tidyverse)


###  Import Single-Omics  ###
dnaM_df <- read_csv(
	# "DNAm_meta/c2cp_trimmed_missMethyl_pValues_20210803.csv"
	"analysis_DNAm/c2cp_trimmed_missMethyl_pValues_20210803.csv"
) %>%
	select(pathway, dnamPval = P.DE, dnamFDR = FDR)

rnaSeq_df <- read_csv(
	# "RNASeq_ROSMAP/ROSMAP_only/c2cp_trimmed_fgsea_pValues_20210810.csv"
	"analysis_RNAseq/c2cp_trimmed_fgsea_pValues_20210810.csv"
) %>%
	select(pathway, size, rnaseqPval = pval, rnaseqFDR = padj)

snp_df <- read_csv(
	# "../../LWtest/GWAS_pathways/RESULT/all_bacon_results_wrangled_20210809.csv"
	"analysis_SNP/all_bacon_results_wrangled_20210809.csv"
) %>%
	select(pathway = TERMS, snpPval = pVal, snpFDR = FDR)


###  Join  ###
multiOmicsAD_df <-
	snp_df %>%
	full_join(dnaM_df, by = "pathway") %>%
	full_join(rnaSeq_df, by = "pathway") %>%
	select(pathway, size, everything())

write_csv(
	multiOmicsAD_df,
	file = "AD_pathway_multiomics_results_20210827.csv"
)

# restart



######  MiniMax Top 10  #######################################################
library(pathwayMultiomics)
library(tidyverse)

multiOmicsAD_df <- read_csv(
	"multiomics/MiniMax/AD_pathway_multiomics_results_20210827.csv"
)

MiniMax(
	pValues_df = multiOmicsAD_df %>%
		select(pathway, ends_with("Pval")),
	orderStat = 2,
	method = "parametric",
	drivers_char = c("SNP", "DNAm", "RNAseq")
) %>%
	mutate(MiniMaxFDR = p.adjust(MiniMaxP, method = "fdr")) %>%
	arrange(MiniMaxFDR) %>%
	slice(1:10) %>%
	pull(pathway)


###  Top few gene sets  ###
"PID_PDGFRB_PATHWAY"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6003857/

"WP_CHEMOKINE_SIGNALING_PATHWAY"
# https://www.frontiersin.org/articles/10.3389/fphar.2019.00622/full

"KEGG_HEMATOPOIETIC_CELL_LINEAGE"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7477654/
# This connection seems thin to me. Maybe this one?
# https://www.nature.com/articles/4001913

"PID_TCR_PATHWAY"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3719061/
# Also included in the previous top 10

"WP_REGULATION_OF_TOLLLIKE_RECEPTOR_SIGNALING_PATHWAY"
# https://doi.org/10.3389/fimmu.2019.01000
# Also included in the previous top 10

"KEGG_CHEMOKINE_SIGNALING_PATHWAY"
# https://www.frontiersin.org/articles/10.3389/fphar.2019.00622/full
# Same paper as WP chemokine signaling pathway

"PID_KIT_PATHWAY"
# Full name: Signaling events mediated by Stem cell factor receptor (c-Kit)
# https://www.sciencedirect.com/science/article/pii/S0163725818302304
# Major player in cancers, so the AD connections are recent

"WP_KIT_RECEPTOR_SIGNALING_PATHWAY"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4759323/

"PID_CXCR4_PATHWAY"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6057706/
# https://www.nature.com/articles/s41398-017-0049-7

"REACTOME_TCR_SIGNALING"
# https://www.nature.com/articles/s41586-019-1895-7
# Also included in the previous top 10



######  Wrangle MiniMax Results  ######################################################
miniMaxResults_df <-
  MiniMax(
  	pValues_df = multiOmicsAD_df %>%
  		select(pathway, ends_with("Pval")),
  	orderStat = 2,
  	method = "parametric",
  	drivers_char = c("SNP", "DNAm", "RNAseq")
  ) %>%
	mutate(MiniMaxQ = p.adjust(MiniMaxP, method = "fdr")) %>%
	left_join(
		multiOmicsAD_df %>%
			select(pathway, size, ends_with("FDR")),
		by = "pathway"
	) %>%
	select(
		pathway, size, ends_with("Pval"), ends_with("FDR"), drivers, everything()
	) %>%
	rename(MiniMaxFDR = MiniMaxQ)

write_csv(
	miniMaxResults_df,
	file = "multiomics/MiniMax/AD_pathway_multiomics_results_wDrivers_20210910.csv"
)
