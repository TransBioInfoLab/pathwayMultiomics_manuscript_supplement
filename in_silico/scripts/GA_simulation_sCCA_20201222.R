# Extract sCCA Code from Pucher et al.
# Gabriel Odom
# 20200904
# UPDATED: 20201201/20201222/20210502

# I'm digging through the functions in IamComparison/src/4_sCCA.R
# NOTE: after reading the help files for MultiCCA and MultiCCA.permute, I think
#   there is an error in the original Pucher et al. code: they use MultiCCA() 
#   *before* calling MultiCCA.permute (the PMA:: examples have this reversed).

# UPDATE 20201201: when we originally attempted sCCA, the power on the highest
#   signal simulation replicate was still awful (~10%). Since then, we've added
#   in the GA RNAseq platform, which has much more shared samples (over 70), so
#   we expect the power to be a bit better here. We are attempting sCCA again.

# UPDATE 20201222: the power is still bad, so we are going to run sCCA for only
#   the most extreme simulation design scenario (80% treated features, +50%
#   signal). But this method is 10+ years old, so I'll give it a pass.

# install.packages("PMA")
library(tidyverse)
library(pathwayPCA)
library(PMA)

cluster_PC <- read_gmt(file = "data/cluster_pathway_collection_20201117.gmt")
#


######  Directory Setup  ######################################################
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

simResFolders <- paste0(
  "simResultsRaw/pathwayPCA_sim_20201201/sim", simulationRunLabels_char
)
dataFiles_char <- c(
  "CNV_partition4_delta0.5.RDS",
  "Prot_partition4_delta0.5.RDS", 
  "RNAseq_partition4_delta0.5.RDS"
)



######  Import Data  ##########################################################
# Treated values
designFiles_char <- paste0(
  simResFolders, "/", "indicatorMatrices", simulationRunLabels_char, "_ls.RDS"
)
simDesign_ls <- readRDS(designFiles_char[1])


# Import the Three Data Sets
dataDF_ls <- map(
  .x = dataFiles_char, # strongest signal
  .f = ~{
    readRDS(paste0(simResFolders[1], "/", .x))
  }
)
names(dataDF_ls) <- c("CNV", "Prot", "RNAseq")



######  Wrangle List of Data Sets  ############################################
sharedSamples_char <- Reduce(
  f = dplyr::intersect,
  x = map(dataDF_ls, "Sample")
)
# 71 shared samples

colVars <- function(mat){
  apply(mat, MARGIN = 2, var)
}

dataMat_ls <- 
  dataDF_ls %>% 
  map(slice, match(Sample, sharedSamples_char)) %>% 
  map(select, -Sample) %>% 
  map(as.matrix) %>% 
  # Remove columns with 0 variance
  map(~{
    .x[, colVars(.x) > 10^-8]
  })
# 1710 genes on all 71 samples

treatedSamps_idx <- unique(simDesign_ls$indicatorMatrices$partition1$i)
resp_lgl <- tibble(
  Sample = simDesign_ls$indicatorDimNames[[1]],
  Treated = Sample %in% simDesign_ls$indicatorDimNames[[1]][treatedSamps_idx]
) %>% 
  # the filter(x %in% y) trick left the samples in the wrong order here, but
  #   the correct order above
  slice(match(sharedSamples_char, Sample)) %>% 
  pull(Treated)

rm(dataDF_ls, sharedSamples_char, treatedSamps_idx)



######  PMA sCCA  #############################################################
system.time(
  permute_out <- MultiCCA.permute(
    xlist = dataMat_ls,
    type = "standard",
    nperms = 100,
    trace = FALSE
  )
)
# 10 perms: 34.03 seconds for 3 platforms with 71 samples and 1710 features each
#   100 permutations: 232.91 seconds for same design

# Now that we have sCCA results, we can find the sCCA vectors with best fit to
#   this data set *and the response*.
out <- MultiCCA(
  xlist = dataMat_ls,
  type = "standard",
  penalty = permute_out$bestpenalties,
  ncomponents = 1,
  ws = permute_out$ws.init,
  trace = FALSE
)

out
str(out)
plot(permute_out)

sCCA_res <- unique(
  union_all(
    colnames(dataMat_ls$CNV)[which(out$ws[[1]] != 0)],
    colnames(dataMat_ls$Prot)[which(out$ws[[2]] != 0)],
    colnames(dataMat_ls$RNAseq)[which(out$ws[[3]] != 0)]
  )
)



######  Analyse Results  ######################################################
# Once we have the significant features, perform an over-representation test
#   using Fisher's Exact method for each pathway
testOverrep = function(features, universe, pathway){
  # Input:
  # - features: character vector of features marked as treated by sCCA
  # - universe: union of all features across all data sets
  # - pathway: character vector of features in the chosen pathway
  
  a = intersect(features, pathway)
  b = setdiff(pathway, features)
  c = intersect(setdiff(universe, pathway), features)
  d = intersect(setdiff(universe, pathway), setdiff(universe, features))
  
  ct =  matrix(c(length(a), length(b), length(c), length(d)),
               nrow = 2, byrow = FALSE,
               dimnames = list("Pathway" = c("Yes", "No"), 
                               "Method" = c("Yes", "No")))
  fp = fisher.test(ct, alternative = "greater")$p.value
  
  return(fp)
}

# Apply this to each pathway
PW_p <- sapply(
  cluster_PC$pathways,
  testOverrep,
  features = sCCA_res,
  universe = colnames(dataMat_ls[[1]])
)
names(PW_p) <- cluster_PC$TERMS


# Tabulate Results
sCCAresults_df <- 
  tibble(
    Pathway = names(PW_p),
    FisherP = PW_p
  ) %>% 
  mutate(
    FisherFDR = p.adjust(PW_p, method = "BH"),
    treated = Pathway %in% simDesign_ls$treatedPathways
  )

sCCAresults_df %>% 
  filter(treated) %>% 
  arrange(FisherP)

# Two of the pathways are significant at rawp < 0.05, none for FDR.


###  AUC  ###
library(pROC)
as.numeric(
  roc(
    response = sCCAresults_df$treated,
    predictor = sCCAresults_df$FisherP,
    levels = as.factor(c(TRUE, FALSE)),
    # Case p-values should be lower than control p-values
    direction = "<"
  )$auc
)



######  Wrap up a Function  ###################################################
###  Helper Functions  ###
colVars <- function(mat) {
  apply(mat, MARGIN = 2, var)
}

testOverrep = function(features, universe, pathway) {
  # Input:
  # - features: character vector of features marked as treated by sCCA
  # - universe: union of all features across all data sets
  # - pathway: character vector of features in the chosen pathway
  
  a = intersect(features, pathway)
  b = setdiff(pathway, features)
  c = intersect(setdiff(universe, pathway), features)
  d = intersect(setdiff(universe, pathway), setdiff(universe, features))
  
  ct =  matrix(c(length(a), length(b), length(c), length(d)),
               nrow = 2, byrow = FALSE,
               dimnames = list("Pathway" = c("Yes", "No"), 
                               "Method" = c("Yes", "No")))
  fp = fisher.test(ct, alternative = "greater")$p.value
  
  return(fp)
}


###  Main Function  ###
RunMultiCCA <- function(simRunLabel, PC) {
  # browser()
  
  ###  Data Import  ###
  simResFolder <- paste0(
    "simResultsRaw/pathwayPCA_sim_20201201/sim", simRunLabel
  )
  dataFiles_char <- c(
    "CNV_partition4_delta0.5.RDS",
    "Prot_partition4_delta0.5.RDS", 
    "RNAseq_partition4_delta0.5.RDS"
  )
  
  # Treated values
  simDesign_ls <- readRDS(
    paste0(
      simResFolder, "/", "indicatorMatrices", simRunLabel, "_ls.RDS"
    )
  )
  # Import the Three Data Sets
  dataDF_ls <- map(
    .x = dataFiles_char,
    .f = ~{ readRDS(paste0(simResFolder, "/", .x)) }
  )
  names(dataDF_ls) <- c("CNV", "Prot", "RNAseq")
  
  
  ###  Data Wrangle  ###
  sharedSamples_char <- Reduce(
    f = dplyr::intersect,
    x = map(dataDF_ls, "Sample")
  )
  dataMat_ls <- 
    dataDF_ls %>% 
    map(slice, match(Sample, sharedSamples_char)) %>% 
    map(select, -Sample) %>% 
    map(as.matrix) %>% 
    # Remove columns with 0 variance
    map(~{
      .x[, colVars(.x) > 10^-8]
    })
  
  rm(dataDF_ls, sharedSamples_char)
  
  
  ###  sCCA  ###
  t0 <- Sys.time()
  permute_out <- MultiCCA.permute(
    xlist = dataMat_ls,
    type = "standard",
    nperms = 100,
    trace = FALSE
  )
  out <- MultiCCA(
    xlist = dataMat_ls,
    type = "standard",
    penalty = permute_out$bestpenalties,
    ncomponents = 1,
    ws = permute_out$ws.init,
    trace = FALSE
  )
  t1 <- Sys.time()
  
  
  ###  Compile Results  ###
  sCCA_res <- unique(
    union_all(
      colnames(dataMat_ls$CNV)[which(out$ws[[1]] != 0)],
      colnames(dataMat_ls$Prot)[which(out$ws[[2]] != 0)],
      colnames(dataMat_ls$RNAseq)[which(out$ws[[3]] != 0)]
    )
  )
  
  PW_p <- sapply(
    cluster_PC$pathways,
    testOverrep,
    features = sCCA_res,
    universe = colnames(dataMat_ls[[1]])
  )
  names(PW_p) <- cluster_PC$TERMS
  
  # Tabulate Results
  sCCAresults_df <- 
    tibble(
      Pathway = names(PW_p),
      FisherP = PW_p
    ) %>% 
    mutate(
      FisherFDR = p.adjust(PW_p, method = "BH"),
      treated = Pathway %in% simDesign_ls$treatedPathways
    )
  
  auc_num <- as.numeric(
    pROC::roc(
      response = sCCAresults_df$treated,
      predictor = sCCAresults_df$FisherP,
      levels = as.factor(c(TRUE, FALSE)),
      # Case p-values should be lower than control p-values
      direction = "<"
    )$auc
  )
  
  
  ###  Save  ###
  results_ls <- list(
    Run = simResFolder,
    multiCCA_permute = permute_out,
    multiCCA_results = out,
    enrichGenes = sCCA_res,
    pathwayPvals_df = sCCAresults_df,
    AUC = auc_num,
    computeTime = t1 - t0
  )
  saveRDS(
    object = results_ls,
    file = paste0(
      "./simAnalysis/small_sim_sCCA_20201222/",
      "results_partition4_delta0.5_run", simRunLabel, ".RDS"
    )
  )
    
}

# Test
RunMultiCCA(simRunLabel = "001", PC = cluster_PC)



######  Apply  ################################################################
###  Walk Through Sim Runs  ###
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

system.time(
  walk(
    .x = simulationRunLabels_char[4:100],
    RunMultiCCA,
    PC = cluster_PC
  )
)
# We expect ~4.3 minutes per run.


###  Import / Extract Results  ###
# UPDATE 2021-05-02: I didn't extract the test size from these results.
slug_char <- "simAnalysis/small_sim_sCCA_20201222/results_partition4_delta0.5"
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

res_df <- 
  map_dfr(
    .x = paste0(slug_char, "_run", simulationRunLabels_char, ".RDS"),
    .f = ~{
      
      x_ls <- readRDS(.x)
      testSize_num <- x_ls$pathwayPvals_df %>%
        filter(!treated) %>%
        mutate(T1err = FisherP < 0.05) %>%
        pull(T1err) %>%
        mean()
      
      tibble(
        Run   = x_ls$Run,
        AUC   = x_ls$AUC,
        TSize = testSize_num,
        time  = x_ls$computeTime
      )
      
    }
  ) %>% 
  mutate(Method = "sCCA")

summary(res_df$AUC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1200  0.4139  0.5122  0.5147  0.6233  0.8467 
summary(res_df$TSize)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.08889 0.08889 0.13444 0.20000 0.28889 
mean(as.numeric(res_df$time))
# 3.40 minutes

write_csv(
  x = res_df,
  file = "simAnalysis/small_sim_sCCA_20201222/sCCA_summary_20210502.csv"
)
