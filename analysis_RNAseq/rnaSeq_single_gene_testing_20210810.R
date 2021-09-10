# Single-Gene RNAseq for ROSMAP Data
# Gabriel Odom and Lanyu Zhang
# 2021-08-10



######  Data and Functions  ###################################################

library(tidyverse)
log2Plus1 <- function(x) { log2(x + 1) }


###  Data  ###
system.time(
	rosmapWide_df <- read_csv(
		"./RNASeq_ROSMAP/ROSMAP_only/rnaSeq_20210805.csv",
		col_types = paste0("c", str_dup("d", 41396))
	)
)
covariates_df <- read_csv(
	"./RNASeq_ROSMAP/ROSMAP_only/covariates_20210805.csv"
) %>%
	mutate(
		ENO2  = log2Plus1(ENO2),
		OLIG2 = log2Plus1(OLIG2),
		CD34  = log2Plus1(CD34),
		CD68  = log2Plus1(CD68),
		GFAP  = log2Plus1(GFAP)
	)


###  Wrapper Function  ###
FitSingleGeneLM <- function(ensID, wide_df, covar_df, trans_fun = identity) {
	# browser()

	###  Create Analysis Dataset  ###
	# This only works because we've already confirmed that the two data sets have
	#   the same row order
	data_df <- wide_df[ , ensID, drop = FALSE]
	colnames(data_df) <- "gene"
	data_df <- cbind(data_df, covar_df)


	###  Fit LM and Extract  ###
	mod <- lm(
		formula = trans_fun(gene) ~ ., data = data_df
	)
	modFit_num <- unname(
		summary(mod)$coefficients["braaksc", c("Estimate", "Pr(>|t|)")]
	)


	###  Return  ###
	tibble(
		ensembl_gene_id = ensID,
		stage_estimate = modFit_num[1],
		stage_pVal = modFit_num[2]
	)

}

# Test
FitSingleGeneLM(
	ensID = colnames(rosmapWide_df)[2],
	wide_df = rosmapWide_df,
	covar_df = covariates_df[ , -(1:2)],
	trans_fun = log2Plus1
)



######  Apply  ################################################################
system.time(
	results_df <- map_dfr(
		.x = colnames(rosmapWide_df)[-1],
		.f = FitSingleGeneLM,
		wide_df = rosmapWide_df,
		covar_df = covariates_df[ , -(1:2)]
	)
)
# 11.019 sec for first 1000; ~11 min for 41k genes

results2_df <-
	results_df %>%
	mutate(
		ensembl_gene_id = str_remove(ensembl_gene_id, pattern = "\\..*")
	) %>%
	arrange(stage_pVal) %>%
	mutate(FDR = p.adjust(stage_pVal, method = "fdr")) %>%
	mutate(FDRsignif = FDR < 0.05)


###  Adding Gene Symbols  ###
geneInfo_df <-
	TCGAbiolinks::get.GRCh.bioMart("hg19") %>%
	select(ensembl_gene_id, symbol = external_gene_name) %>%
	unique()

results3_df <-
	results2_df %>%
	left_join(geneInfo_df, by = "ensembl_gene_id") %>%
	select(ensembl_gene_id, symbol, everything())


write_csv(
	results3_df,
	file = "./RNASeq_ROSMAP/ROSMAP_only/rnaSeq_single_gene_pVals_20210810.csv"
)

# Restart



######  fgsea Analysis  #######################################################

###  Setup and Data  ###
library(fgsea)
library(biomaRt)
library(tidyverse)

# results
results_df <- read_csv(
	file = "./RNASeq_ROSMAP/ROSMAP_only/rnaSeq_single_gene_pVals_20210810.csv"
)

# map from EnsemblID to EntrezID
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes_df <- getBM(
	filters = "ensembl_gene_id",
	attributes = c("ensembl_gene_id", "entrezgene_id"),
	values = results_df$ensembl_gene_id,
	mart = mart
)
# We match 38470 of the 41396 targets
results2_df <-
	genes_df %>%
	rename(entrez_gene_id = entrezgene_id) %>%
	left_join(
		results_df %>%
			select(ensembl_gene_id, symbol, stage_pVal),
		by = "ensembl_gene_id"
	) %>%
	mutate(entrez_gene_id = as.character(entrez_gene_id)) %>%
	as_tibble()

rm(mart, genes_df, results_df)


# our PC is in EntrezIDs, not Ensembl IDs.
results_num <- -log10(results2_df$stage_pVal)
names(results_num) <- results2_df$entrez_gene_id


# gene sets
c2cpTrim_PC <- pathwayPCA::read_gmt(
	file = "../../gene_sets/c2_cp_entrez_trimmed_20210709.gmt"
)

# fgsea requires the pathway collection be in a slightly different format
c2cpTrim_ls <- c2cpTrim_PC$pathways
names(c2cpTrim_ls) <- c2cpTrim_PC$TERMS
rm(c2cpTrim_PC)


###  Gene Set Analysis  ###
set.seed(12345)
system.time(
	res_fgsea <- fgsea(
		pathways = c2cpTrim_ls,
		stats = results_num,
		scoreType = "pos"
	)
)
# 4.776 sec for 38470 genes in 2833 pathways

write_csv(
	res_fgsea,
	file = "./RNASeq_ROSMAP/ROSMAP_only/c2cp_trimmed_fgsea_pValues_20210810.csv"
)


