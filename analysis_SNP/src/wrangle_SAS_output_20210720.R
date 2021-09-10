# Wrangle Results
# Gabriel Odom
# 2021-07-20



######  Overview and Setup  ###################################################
# Lily finished running the SAS code to do the pathway-level SNP analysis. We
#   need to import this results file, calculate the pathway p-values (including
#   for pathways that the algorithm did not converge), match the pathway names
#   (SAS truncated the pathway names to 50 characters), and save the cleaned
#   results

library(readxl)
library(pathwayPCA)
library(tidyverse)

sasOutRaw_df <- read_csv(
	"../LWtest/GWAS_pathways/RESULT/all_results_bacon_correction_7-21-2021.csv"
)

c2cpTrim_PC <- read_gmt(
	file = "../gene_sets/c2_cp_entrez_trimmed_20210709.gmt",
	description = TRUE
)



######  Pathway Names  ########################################################

sasOutRaw_df <-
	sasOutRaw_df %>%
	arrange(geneset) %>%
	mutate(TERMS = c2cpTrim_PC$TERMS)

sasOutClean_df <-
	sasOutRaw_df %>%
	mutate(converged = str_detect(Reason, pattern = "satisfied")) %>%
	select(
		TERMS, converged, tValue, DF, pVal = pValue.bacon, FDR = fdr.bacon
	) %>%
	mutate(
		pVal = case_when(
			converged ~ pVal,
			!converged ~ 1
		)
	) %>%
	select(-converged) %>%
	ungroup() %>%
	mutate(pAdjust = p.adjust(pVal, method = "bonferroni"))



######  Add Collections  ######################################################
sasOut2_df <-
	sasOutClean_df %>%
	rowwise() %>%
	mutate(
		collection = str_split(TERMS, pattern = "_", simplify = TRUE)[, 1]
	) %>%
	mutate(signif = pAdjust < 0.05) %>%
	ungroup()

write_csv(
	sasOut2_df,
	file = "../LWtest/GWAS_pathways/RESULT/all_bacon_results_wrangled_20210809.csv"
)
