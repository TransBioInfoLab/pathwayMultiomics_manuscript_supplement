# Cluster Shared Genes into Independent Pathways
# Gabriel Odom
# 2020-11-17

library(tidyverse)
library(pathwayPCA)



######  Colon Cancer Cases  ###################################################
# Grab data from the TCGA COADREAD page on:
# http://www.linkedomics.org/data_download/TCGA-COADREAD/


# RNAseq (GenomaAnalyzer)
coadRNAseqGA_df <- read_delim(
  "data/Human__TCGA_COADREAD__UNC__RNAseq__GA_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadRNAseqGA_df <- TransposeAssay(coadRNAseqGA_df)
# 222 samples, 6149 features


# Copy Number Variation (Log-Ratio)
coadCNV_df <- read_delim(
  "data/Human__TCGA_COADREAD__BI__SCNA__SNP_6.0__01_28_2016__BI__Gene__Firehose_GISTIC2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadCNV_df <- TransposeAssay(coadCNV_df)
# 616 samples, 24776 features


# Proteomics
coadProt_df <- read_delim(
  "data/Human__TCGA_COADREAD__VU__Proteome__Velos__01_28_2016__VU__Gene__CDAP_UnsharedPrecursorArea_r2.cct", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadProt_df <- TransposeAssay(coadProt_df)
# 90 samples, 5538 features


# Clinical
coadClinical_df <- read_delim(
  "data/Human__TCGA_COADREAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadClinical_df <-
  TransposeAssay(coadClinical_df) %>% 
  select(Sample, overall_survival, status) %>% 
  drop_na()


###  Check Sample Overlaps  ###
# SEE scripts/GA_simulation_pt1_20201117.R


######  Compose Shared Gene Dataset  ##########################################
commonGenes_char <- readRDS(file = "data/common_genes_20201117.RDS") # 1710
matchGene_ls <- list(
  cnv  = coadCNV_df %>% 
    select(one_of(commonGenes_char)),
  rna  = coadRNAseqGA_df %>% 
    select(one_of(commonGenes_char)),
  prot = coadProt_df %>% 
    select(one_of(commonGenes_char))
)

coadMatchedGenes_df <- 
  matchGene_ls %>% 
  map(
    .f = ~{
      mutate_all(.x, function(x) {as.numeric(scale(x))} )
    }
  ) %>% 
  bind_rows()
# 928 samples, 1710 features

rm(matchGene_ls)



######  Clustering  ###########################################################
# We want to find out how far apart the genes are, so we must transpose our
#   data:
coadMatchedGenesT_mat <- t(coadMatchedGenes_df)

system.time(
  allGenesEuc_dist <- dist(coadMatchedGenesT_mat, method = "euclidean")
)
# Euclidean: 12.147 seconds (Mac)

geneClust <- hclust(d = allGenesEuc_dist, method = "ward.D")

plot(geneClust)

geneClust %>% 
  cutree(k = 65) %>% 
  table() %>% 
  as.matrix() %>% 
  as.vector() %>% 
  summary()

# Using EUCLIDEAN + ward.D (k = 200):
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    5.00    8.00    8.55   11.00   27.00 
# Using EUCLIDEAN + WARD (k = 100):
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3.0    11.0    17.0    17.1    22.0    47.0
# Using EUCLIDEAN + WARD (k = 50):
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   9.0    22.0    33.0    34.2    44.5    74.0

# Now let's extract the pathways:
geneClust
groupMap_df <- as_tibble(
  cbind(
    gene  = rownames(coadMatchedGenesT_mat),
    group = cutree(geneClust, k = 50)
  )
)

pathways_ls <- map(
  .x = unique(groupMap_df$group),
  .f = ~{
    groupMap_df %>% 
      filter(group == .x) %>% 
      pull(gene)
  }
)

cluster_PC <- CreatePathwayCollection(
  sets_ls = pathways_ls,
  TERMS = paste0(
    "cluster",
    str_pad(seq_along(pathways_ls), width = 2, pad = "0")
  )
)
# 50 pathways that cover 1710 genes with median cardinality of 33 genes and
#   range [9, 74]
write_gmt(cluster_PC, file = "data/cluster_pathway_collection_20201117.gmt")
