# MiniMax Statistic Under the Null
# Gabriel Odom
# 2020-12-02



#######  Introduction  ########################################################
# We want to find the MiniMax statistic for each pathway under the null
#   distribution (i.e., no signal added to the original data, only the clinical
#   response has been randomly assigned). These statistic values will give us
#   the information we need to later estimate the distribution of the MiniMax
#   statistic and calculate its test size.   

library(tidyverse)
library(pathwayPCA)
commonGenes_char <- readRDS(file = "data/common_genes_20201117.RDS")
cluster_PC <- read_gmt(file = "data/cluster_pathway_collection_20201117.gmt")



######  Import + Wrangle Data  ################################################
# Grab data from the TCGA COADREAD page on:
# http://www.linkedomics.org/data_download/TCGA-COADREAD/



###  RNAseq (GA)  ###
coadRNAseqGA_df <- read_delim(
  "data/Human__TCGA_COADREAD__UNC__RNAseq__GA_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)

# There are missing values in this data set (139 of 1,365,078 values):
summary(as.vector(as.matrix(coadRNAseqGA_df[, -1])))
coadRNAseqGAimpute_ls <- impute::impute.knn(as.matrix(coadRNAseqGA_df[, -1]))
summary(as.vector(coadRNAseqGAimpute_ls$data))
# We've imputed these missing values with default parameters to impute.knn()

coadRNAseqGAimpute_df <-
  data.frame(
    attrib_name = coadRNAseqGA_df[, 1],
    as.data.frame(coadRNAseqGAimpute_ls$data)
  ) %>% 
  TransposeAssay() %>% 
  select(Sample, one_of(commonGenes_char)) %>% 
  mutate_if(is.numeric, scale)

rm(coadRNAseqGA_df, coadRNAseqGAimpute_ls)


###  Copy Number Variation (Log-Ratio)  ###
coadCNV_df <- read_delim(
  "data/Human__TCGA_COADREAD__BI__SCNA__SNP_6.0__01_28_2016__BI__Gene__Firehose_GISTIC2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
) %>% 
  filter(attrib_name %in% commonGenes_char) %>% 
  TransposeAssay() %>% 
  mutate_if(is.numeric, scale)


###  Proteomics  ###
coadProt_df <- read_delim(
  "data/Human__TCGA_COADREAD__VU__Proteome__Velos__01_28_2016__VU__Gene__CDAP_UnsharedPrecursorArea_r2.cct", 
  "\t", escape_double = FALSE, trim_ws = TRUE
) %>% 
  filter(attrib_name %in% commonGenes_char) %>% 
  TransposeAssay() %>% 
  mutate_if(is.numeric, scale)


###  Clinical  ###
coadClinical_df <- read_delim(
  "data/Human__TCGA_COADREAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "\t", escape_double = FALSE, trim_ws = TRUE
) %>% 
  TransposeAssay() %>% 
  select(Sample, overall_survival, status) %>% 
  drop_na()



######  PathwayPCA Wrapper  ###################################################
results_ls <- vector(mode = "list", length = 100L)

t0 <- Sys.time()
results_ls[26:100] <- map(
  .x = 26:100,
  .f = ~{
    
    ###  Omics* Containers  ###
    coadClinicalRand_df <- tibble(
      Sample = coadClinical_df$Sample,
      survTime = sample(coadClinical_df$overall_survival)
    )
    
    CNV_omics <- CreateOmics(
      assayData_df = coadCNV_df,
      pathwayCollection_ls = cluster_PC,
      response = coadClinicalRand_df,
      respType = "regression"
    )
    
    RNAseq_omics <- CreateOmics(
      assayData_df = coadRNAseqGAimpute_df,
      pathwayCollection_ls = cluster_PC,
      response = coadClinicalRand_df,
      respType = "regression"
    )
    
    prot_omics <- CreateOmics(
      assayData_df = coadProt_df,
      pathwayCollection_ls = cluster_PC,
      response = coadClinicalRand_df,
      respType = "regression"
    )
    
    
    ###  AES-PCA Calls  ###
    CNV_aespcOut <- AESPCA_pVals(
      object = CNV_omics,
      parallel = TRUE,
      numCores = 12,
      adjustment = "BH"
    )
    RNAseq_aespcOut <- AESPCA_pVals(
      object = RNAseq_omics,
      parallel = TRUE,
      numCores = 12,
      adjustment = "BH"
    )
    prot_aespcOut <- AESPCA_pVals(
      object = prot_omics,
      parallel = TRUE,
      numCores = 12,
      adjustment = "BH"
    )
    
    
    ###  Extract p-Values  ###
    allPvals_df <- 
      left_join(
        CNV_aespcOut$pVals_df %>% 
          select(pathways, rawp) %>% 
          rename(CNV_p = rawp),
        RNAseq_aespcOut$pVals_df %>% 
          select(pathways, rawp) %>% 
          rename(RNA_p = rawp)
      ) %>% 
      left_join(
        prot_aespcOut$pVals_df %>%
          select(pathways, rawp) %>% 
          rename(Prot_p = rawp)
      ) %>% 
      rowwise() %>% 
      mutate(
        MiniMax = sort(
          c_across(where(is.numeric)),
          decreasing = FALSE
        )[2]
      ) %>% 
      arrange(MiniMax) %>% 
      ungroup()
    
    
    ###  Return  ###
    print(paste("Replicate", .x, "complete at", Sys.time()))
    allPvals_df
    
  }
)
t1 <- Sys.time()
# 10.92502 min for first 5 replicates; 42.30286 min for next 20; 2.580547 hrs
#   for last 75
# results2_ls <- map(results_ls, ungroup)

saveRDS(
  results_ls,
  file = "simResultsRaw/pathwayPCA_sim_20201201/miniMax_test_size_20201201.RDS"
)



######  Inspect Results  ######################################################
map(results_ls, "CNV_p") %>% unlist %>% hist()
map(results_ls, "RNA_p") %>% unlist %>% hist()
map(results_ls, "Prot_p") %>% unlist %>% hist()


