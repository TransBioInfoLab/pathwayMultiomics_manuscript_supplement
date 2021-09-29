# Create Simulated Data Sets
# Gabriel Odom
# 2020-08-26


library(tidyverse)
library(pathwayPCA)
library(slam)



######  Base Simulation Layer  ################################################
# In the previous step (detail shown in scripts/GA_simulation_pt2_20201201.R),
#   we generated 100 folders for the 100 simulation runs. In those folders, we
#   also saved 5 indicator matrices of dimension 631 x 1710 (631 total samples
#   from TCGA COADREAD data and 1710 features shared across proteomics, RNAseq
#   from the GenomeAnalyzer platform, and CNV.)
# Now, we will mark the location of these indicator matrices, and import the
#   base genomic data that will be combined with the indicator matrices to make
#   the treated data.
# Grab data from the TCGA COADREAD page on:
# http://www.linkedomics.org/data_download/TCGA-COADREAD/


###  Folders  ###
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

simResFolders <- paste0(
  "simResultsRaw/pathwayPCA_sim_20201201/sim", simulationRunLabels_char
)


###  Genes and Indicators  ###
commonGenes_char <- readRDS(file = "data/common_genes_20201117.RDS")
simResFiles_char <- paste0(
  simResFolders,
  "/indicatorMatrices",
  simulationRunLabels_char,
  "_ls.RDS"
)


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
)
coadCNV_df <-
  TransposeAssay(coadCNV_df) %>% 
  select(Sample, one_of(commonGenes_char)) %>% 
  mutate_if(is.numeric, scale)


###  Proteomics  ###
coadProt_df <- read_delim(
  "data/Human__TCGA_COADREAD__VU__Proteome__Velos__01_28_2016__VU__Gene__CDAP_UnsharedPrecursorArea_r2.cct", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadProt_df <-
  TransposeAssay(coadProt_df) %>% 
  select(Sample, one_of(commonGenes_char)) %>% 
  mutate_if(is.numeric, scale)


###  Clinical  ###
coadClinical_df <- read_delim(
  "data/Human__TCGA_COADREAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadClinical_df <-
  TransposeAssay(coadClinical_df) %>% 
  select(Sample, overall_survival, status) %>% 
  drop_na()



######  Data Generation Function  #############################################
# In each of these folders, we have 4 indicator matrices. We need a function
#   that will take in the three real data sets, create the 3 x 4 x 5 simulation
#   data sets, then save these data sets in their appropriate folders.

SaveSimData <- function(indicatorPath_char,
                        dir_char,
                        df, dfType_char,
                        effectSize = c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  # browser()
  
  indicator_ls <- readRDS(file = indicatorPath_char)
  keepRows_idx <- match(df[["Sample"]], indicator_ls$indicatorDimNames[[1]])
  orderedGenes_char <- indicator_ls$indicatorDimNames[[2]]
  
  # If the columns aren't in the same order, adding two matrices is meaningless
  df2 <- df[, c("Sample", orderedGenes_char)]
  
  indicatorMat_ls <-
    map(
      .x = indicator_ls$indicatorMatrices,
      .f = ~{
        
        idx_mat <- as.matrix(.x)
        dimnames(idx_mat) <- indicator_ls$indicatorDimNames
        idx_mat[keepRows_idx , ]
        
      }
    ) 
  
  designCombos_char <- as.vector(
    outer(
      names(indicator_ls$indicatorMatrices),
      effectSize,
      paste, sep = "_"
    )
  )
  designs_ls <- str_split(designCombos_char, pattern = "_")
  names(designs_ls) <- designCombos_char
  
  walk(
    .x = designs_ls,
    .f = ~{
      
      # browser()
      
      out_mat <- 
        as.matrix(df2[, -1]) + 
          indicatorMat_ls[[.x[1]]] * as.numeric(.x[2])
      
      saveRDS(
        object = data.frame(
          Sample = df[["Sample"]],
          out_mat,
          stringsAsFactors = FALSE
        ),
        file = paste0(
          dir_char, "/", dfType_char, "_", .x[[1]], "_delta", .x[[2]], ".RDS"
        )
      )
      
    }
  )
  
}


# Test
SaveSimData(
  indicatorPath_char = simResFiles_char[1],
  dir_char = simResFolders[1],
  df = coadRNAseqGAimpute_df,
  dfType_char = "RNAseq"
)

table(
  round(
    as.vector(
      as.matrix(RNAseq_partition4_delta0.1[,-1]) - 
        as.matrix(
          coadRNAseqGAimpute_df[, colnames(RNAseq_partition4_delta0.1[,-1])]
        )
    ),
    digits = 12
  )
)
# There are a few (about 40) entries that differ at machine
#   precision (4.024558e-16 or less). Therefore, I've added the round() call.


###  Apply  ###
system.time(
  walk(
    .x = seq_len(100),
    .f = ~{
      SaveSimData(
        indicatorPath_char = simResFiles_char[.x],
        dir_char = simResFolders[.x],
        # Update the data types here:
        df = coadRNAseqGAimpute_df, dfType_char = "RNAseq"
      )
    }
  )
)

# 12.6 min for CNV on FIU PC; 7.8 min for RNAseq on FIU PC; 3.4 min for Prot

