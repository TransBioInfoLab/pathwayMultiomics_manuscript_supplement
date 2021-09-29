# Attempting Multi-Omics Simulation Data from TCGA Cases
# Gabriel Odom
# 2020-12-01



######  Overview  #############################################################
# We're taking TCGA cases (from colon cancer) to simulate treated and untreated
#   samples. To set up the simulation study, we will
#   1. Identify 50 pathways with no overlap, p_1, ..., p_50. We should create
#      pathways with at least 5 genes in them (20% of 5 is 1; we can show
#      treatment proportion as low as 20%). 
#      NOTE: this is complete: "cluster_pathway_collection_20201117.gmt" 
#   2. Identify cases (not necessarily matched) with CNV, protein, and RNAseq
#      values, n_CNV, n_prot, n_gene. Using cases will ensure that the
#      correlation structures will be the same for all data. Our technique
#      doesn't necessarily pick out second-moment differences.
#      NOTE: We will use the COAD data set with the genes listed in 
#      "common_genes_20201117.RDS"
#   3. Standardise all features to have mean 0 and standard deviation 1.
#   4. Ensure we have complete survival information for each subject (if we 
#      potentially would like to analyse survival response in addition to binary
#      classification).

# Given this pathway collection and these three data sets, we will simulate a
#   "subtype" of this cancer by treating half of the samples at random (ensuring
#   to treat those samples which are matched as matched. In order to do this, 
#   we should find the samples that are matched across pairs of the three data
#   sets and treat half of the samples in the intersection, then treat samples
#   outside these intersections until half of the samples from each data set
#   have been treated. This will require a table of treated and untreated sample
#   IDs for the three platforms for each simulation run.) In order to simulate
#   these data sets, we should use the simDataRaw/ directory. We can create
#   subdirectories by design point and sub-subdirectories by run ID. This will
#   have similar design to the "IamComparison" project of Pucher et al. (2019).

# Concerning the design points, we would like control over the following:
#   1. How many pathways of the 50 will be disregulated for the cancer subtype?
#      We will default to 5.
#   2. What proportion of genes in each disregulated pathway will be treated? We
#      will have this range over 20%, 40%, ..., 80%. 
#   3. What will the treatment effect be? Use delta = +0.1x, +0.2x, ..., +0.5x
#      of the standard deviation.


library(tidyverse)
library(pathwayPCA)



######  Colon Cancer Cases  ###################################################
# Grab data from the TCGA COADREAD page on:
# http://www.linkedomics.org/data_download/TCGA-COADREAD/

commonGenes_char <- readRDS(file = "data/common_genes_20201117.RDS")


###  Import and Transpose  ###
# RNAseq (GA)
coadRNAseqGA_df <- read_delim(
  "data/Human__TCGA_COADREAD__UNC__RNAseq__GA_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadRNAseqGA_df <-
  TransposeAssay(coadRNAseqGA_df) %>% 
  select(Sample, one_of(commonGenes_char)) %>% 
  mutate_if(is.numeric, scale)


# Copy Number Variation (Log-Ratio)
coadCNV_df <- read_delim(
  "data/Human__TCGA_COADREAD__BI__SCNA__SNP_6.0__01_28_2016__BI__Gene__Firehose_GISTIC2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadCNV_df <-
  TransposeAssay(coadCNV_df) %>% 
  select(Sample, one_of(commonGenes_char)) %>% 
  mutate_if(is.numeric, scale)


# Proteomics
coadProt_df <- read_delim(
  "data/Human__TCGA_COADREAD__VU__Proteome__Velos__01_28_2016__VU__Gene__CDAP_UnsharedPrecursorArea_r2.cct", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadProt_df <-
  TransposeAssay(coadProt_df) %>% 
  select(Sample, one_of(commonGenes_char)) %>% 
  mutate_if(is.numeric, scale)


# Clinical
coadClinical_df <- read_delim(
  "data/Human__TCGA_COADREAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadClinical_df <-
  TransposeAssay(coadClinical_df) %>% 
  select(Sample, overall_survival, status) %>% 
  drop_na()


###  Table of Sample Overlaps  ###
allSamples_char <- 
  union_all(
    coadClinical_df$Sample,
    coadCNV_df$Sample,
    coadRNAseqGA_df$Sample,
    coadProt_df$Sample
  ) %>% 
  unique()

# See Venn diagram in "GA_simulation_pt1_20201117.R"



######  Create Design Template  ###############################################
# NOTE 20201201: because the overall union of samples is the same (see the
#   "Create Design Template" section of "GA_simulation_pt1_20201117.R"), we can
#   skip this work (but I include the code for records).

###  Subtypes by Run  ###
# simulatedSubtypes_char <- rep(c("A", "B"), length = nrow(allSamples_df))
# simulationRunLabels_char <- 
#   seq_len(100) %>% 
#   as.character() %>% 
#   str_pad(width = 3, pad = "0")
# 
# set.seed(12345)
# simDesignTest_df <- map_dfr(
#   .x = simulationRunLabels_char,
#   .f = ~{
#     tibble(
#       Sample   = allSamples_df$Sample,
#       simRun   = .x,
#       simType  = sample(simulatedSubtypes_char, replace = FALSE)
#     )
#   }
# )
# write_csv(
#   simDesignTest_df,
#   path = "simDataRaw/treated_samples_map_20200720.csv"
# )

simDesignTest_df <- read_csv(
  file = "simDataRaw/treated_samples_map_20200720.csv"
)

# What are the proportions for proteomics? (The other two platforms have much
#   larger sample sizes, so I'm not worried about them being slightly unbalanced)
simDesignTest_df %>% 
  filter(Sample %in% coadProt_df$Sample) %>% 
  group_by(simRun) %>% 
  summarise(propA = mean(simType == "A")) %>% 
  pull(propA) %>% 
  quantile(seq(0, 1, 0.1))
# That looks pretty good to me: 80% of the data between 0.43 and 0.57. 100% of
#   the sample balances are between 0.37 and 0.69.

# Just to cover my bases, the other two platforms are
simDesignTest_df %>% 
  filter(Sample %in% coadRNAseqGA_df$Sample) %>% 
  group_by(simRun) %>% 
  summarise(propA = mean(simType == "A")) %>% 
  pull(propA) %>% 
  quantile(seq(0, 1, 0.1))
# 80% between 0.46 and 0.54; 100% between 0.42 and 0.59
simDesignTest_df %>% 
  filter(Sample %in% coadCNV_df$Sample) %>% 
  group_by(simRun) %>% 
  summarise(propA = mean(simType == "A")) %>% 
  pull(propA) %>% 
  quantile(seq(0, 1, 0.1))
# Perfect, just perfect.

rm(list = ls())



######  Simulation Details  ###################################################
###  Design Points by Run  ###
# Now we will make a data frame of all combinations of design points. These are
#   - the proportion of genes in the pathway which we will treat: 20%, ..., 80%
#     (by 20% steps)
#   - the signal strength in terms of standard deviation: 0.1x, 0.2x, ..., 0.5x.
#     Because the features are all standardised, we will add 0.1, 0.2, ..., 0.5.
# This will result in 4 x 5 x 100 data sets. I need to save the simulation run
#   (which will tell me which samples get treated), which pathways will contain
#   the treated genes, the proportion of treated genes for that run, which
#   genes will be treated, and the strength of that treatment. In order to
#   ensure an "apples to apples" comparison, we need the genes selected at a 
#   smaller treatment proportion to be a proper subset of the genes selected at
#   a higher proportion.

###  Steps  ###
# For each simulation run:
#   a) mark the samples in the treatment group
#   b) mark the 5 pathways to be treated
#   c) for each pathway:
#      - partition the pathway into 5 parts
#      - label the partitions 1-5
#      - for 80% treated genes, take the genes in partitions 1-4; for 60%, take
#        the genes in partitions 1-3; and so on until for 20% treated genes, we
#        take partition 1. This means we only have to store the genes once per
#        pathway per simulation run.
#   d) for each set of partitions (1-4, 1-3, 1-2, and 1):
#      - create an indicator matrix (1 for treatment group and treated gene, 0
#        otherwise)
#      - save these indicator matrices in sparse matrix form
#   e) save a list of: a character string of the 5 pathway names; four sparse
#      logical matrices; the (shared) row and column names of these matrices.
#      This format will allow us to add or delete treated pathways or treatment
#      subjects while preserving an "apples to apples" comparison set.

# Once we have these 100 lists, we create simulation data by the following. For
#   each platform p = 1, 2, 3, each simulation replicate i = 1, ..., 100, for
#   each treatment proportion j = 1, 2, 3, 4, row-match each treatment
#   indicator matrix I_{i,j} to the data set of choice X_p, yielding I_{p,i,j}.
#   Then, for each treatment level d = 1, 2, 3, 4, 5, define the simulated
#   data to be X_{p,i,j,d} := X_p + I_{p,i,j} * d.

""


######  Function Guts  ########################################################
simDesignTest_df <- read_csv("simDataRaw/treated_samples_map_20200720.csv")
cluster_PC <- read_gmt(file = "data/cluster_pathway_collection_20201117.gmt")

library(slam)

IndicateTreatment <- function(simResp_df, PC,
                              nTreatedPaths = 5L,
                              nPartitions = 5L,
                              maxTreatPart = 4L,
                              path = NULL){
  # browser()
  
  ###  Set the Treated Genes  ###
  # no reason to do this inside the function other than keeping the global
  #   environment clean
  allGenes_char <- unique(
    unlist(PC$pathways)
  )
  
  # Randomly assign pathways to be treated
  treatedPaths_idx <- sample(length(PC$TERMS), size = nTreatedPaths)
  treatedPaths_ls <- PC$pathways[treatedPaths_idx]
  names(treatedPaths_ls) <- PC$TERMS[treatedPaths_idx]
  
  geneTreatmentGroups_df <- do.call(
    rbind,
    lapply(
      X = treatedPaths_ls,
      FUN = function(path){
        data.frame(
          # Some pathways have all the genes in the same family in series
          gene = sample(path),
          treatLVL = rep_len(
            seq_len(nPartitions),
            length.out = length(path)
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  # The rownames are the cluster name with ".NN" appended to them (where "NN"
  #   is some integer); e.g., "cluster47.21". We then split these strings into
  #   two columns of a character matrix and keep the first column as an atomic
  #   character vector.
  geneTreatmentGroups_df$pathway <- 
    str_split_fixed(rownames(geneTreatmentGroups_df), "\\.", n = 2)[, 1]
  rownames(geneTreatmentGroups_df) <- NULL
  
  
  ###  Indicator Matrices  ###
  # Create an indicator matrix for samples (given a treated sample, mark all
  #   genes with a 1; depending on the design, this should be a matrix half
  #   full of 1s)
  treatedSamples_mat <- matrix(
    data = simResp_df$simType == "B",
    nrow = nrow(simResp_df),
    ncol = length(allGenes_char),
    byrow = FALSE
  )
  
  # Create an indicator matrix for genes at each treatment level (given a 
  #   treated gene, mark all samples with a 1; depending on the design, these 5
  #   matrices should have between 1% - 10% 1s)
  treatedGenesMat_ls <- lapply(
    X = seq_len(maxTreatPart),
    FUN = function(i) {
      
      whichGeneRows_idx <- geneTreatmentGroups_df$treatLVL %in% seq_len(i)
      whichGenes_char <- geneTreatmentGroups_df$gene[whichGeneRows_idx]
      
      matrix(
        data = allGenes_char %in% whichGenes_char + 0,
        nrow = nrow(simResp_df),
        ncol = length(allGenes_char),
        byrow = TRUE
      )
    }
  )
  names(treatedGenesMat_ls) <- paste0("partition", seq_len(maxTreatPart))
  
  # Multiply
  indicatorMatrices_ls <- lapply(
    X = treatedGenesMat_ls,
    FUN = function(x) { 
      as.simple_triplet_matrix(x * treatedSamples_mat)
    }
  )
  
  
  ###  Return  ###
  data_ls <- list(
    treatedPathways = names(treatedPaths_ls),
    indicatorDimNames = list(simResp_df$Sample, allGenes_char),
    indicatorMatrices = indicatorMatrices_ls
  )
  
  if (is.null(path)) {
    data_ls
  } else {
    saveRDS(object = data_ls, file = path)
  }
  
}

# Test
test_ls <- IndicateTreatment(
  simResp_df = simDesignTest_df %>% 
    filter(simRun == "001"),
  PC = cluster_PC
)
rm(test_ls)

IndicateTreatment(
  simResp_df = simDesignTest_df %>% 
    filter(simRun == "001"),
  PC = cluster_PC,
  path = paste0("test", "001", ".RDS")
)



######  Apply  ################################################################
# We will create these files on the FIU PC. I don't think I'll have space on 
#   the UM PC, and I know I don't have space on my Mac.
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

simResFolders <- paste0(
  "simResultsRaw/pathwayPCA_sim_20201201/sim", simulationRunLabels_char
)
walk(.x = simResFolders, .f = dir.create)

walk2(
  .x = simulationRunLabels_char,
  .y = simResFolders,
  .f = ~{
    IndicateTreatment(
      simResp_df = simDesignTest_df %>% 
        filter(simRun == .x),
      PC = cluster_PC,
      path = paste0(
        .y, paste0("/indicatorMatrices", .x, "_ls.RDS")
      )
    )
  }
)


