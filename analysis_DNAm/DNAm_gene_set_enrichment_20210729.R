# Apply missMethyl to AD Data
# Gabriel Odom
# 2021-07-29



######  Extract Single Probes and Combine  ####################################
library(coMethDMR)

StackRegionCpGs <- function(regionName, genome = "hg19", arrayType = "450k"){

	regionCpGs_char <- coMethDMR::GetCpGsInRegion(
		regionName_char = regionName,
		genome = genome, arrayType = arrayType
	)

	data.frame(
		inputRegion = rep_len(regionName, length.out = length(regionCpGs_char)),
		cpg = regionCpGs_char,
		stringsAsFactors = FALSE
	)

}

# Data set
signifDMRs_df <- read.csv(
	file = "./DNAm_meta/meta_analysis_no_crossHyb_smoking_ov_comb_p.csv"
)

# I don't want to load the Tidyverse because I'm already using Bioconductor
#   packages; sometimes those two ecosystems don't play nice.
system.time(
	regionSignifCpGs_ls <- lapply(signifDMRs_df$inputRegion, StackRegionCpGs)
)
# 24.981 sec for 119 regions
regionSignifCpGs_df <- do.call(rbind, regionSignifCpGs_ls)
# 742 probes within regions

write.csv(
	merge(
		x = regionSignifCpGs_df, y = signifDMRs_df,
		by = "inputRegion"
	),
	file = "./DNAm_meta/meta_analysis_no_crossHyb_smoking_ov_comb_p_singleCpG_20210803.csv",
	row.names = FALSE
)

# Restart here



######  Gene Set Enrichment  ##################################################

###  Setup  ###
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# As of 20210803, the package above caused missMethyl to fail in loading, so I'm
#   loading it separately
library(missMethyl)


###  Data and Significant CpGs  ###
individSignifCpGs_df <- read.csv(
	file = "./DNAm_meta/meta_analysis_single_cpg_sig_no_crossHyb_smoking_with_state_greatAnnot_df.csv"
)
regionSignifCpGs_df <- read.csv(
	file = "./DNAm_meta/meta_analysis_no_crossHyb_smoking_ov_comb_p_singleCpG_20210803.csv"
)
signifCpGs_char <- c(individSignifCpGs_df$cpg, regionSignifCpGs_df$cpg)


###  Pathway Collection  ###
c2cpTrim_PC <- pathwayPCA::read_gmt(
	file = "../../gene_sets/c2_cp_entrez_trimmed_20210709.gmt"
)

# missMethyl requires the pathway collection be in a slightly different format
c2cpTrim_ls <- c2cpTrim_PC$pathways
names(c2cpTrim_ls) <- c2cpTrim_PC$TERMS

rm(c2cpTrim_PC)


###  missMethyl  ###
system.time(
	c2cpTrimEnrichment_df <- gsameth(
		sig.cpg = signifCpGs_char,
		collection = c2cpTrim_ls,
		array.type = "450K"
	)
)
# 17.749 sec for 2833 C2CP

c2cpTrimEnrichment2_df <- data.frame(
	pathway = rownames(c2cpTrimEnrichment_df),
	c2cpTrimEnrichment_df
)
rownames(c2cpTrimEnrichment2_df) <- NULL

write.csv(
	c2cpTrimEnrichment2_df,
	file = "./DNAm_meta/c2cp_trimmed_missMethyl_pValues_20210803.csv",
	row.names = FALSE
)
