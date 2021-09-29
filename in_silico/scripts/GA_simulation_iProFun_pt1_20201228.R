# Analyse Simulated Data with iProFun::
# Gabriel Odom
# 2020-12-28

# Lily reminded me of Pei's package for integrative multi-omics. I'm going to
#   test it now. They do not have a vignette, but they do have a README:
# https://github.com/songxiaoyu/iProFun

# UPDATE 20201228: the power isn't bad for the strongest signal design point
#   (Median = 0.7733; Q3 = 0.8517; Max = 0.9667). 

# install.packages("metRology")
# devtools::install_github("songxiaoyu/iProFun")

library(tidyverse)
library(pathwayPCA)
library(iProFun)

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
dataFiles_char <- list.files(simResFolders[2])

simAnalysisFolders <- paste0(
  "simAnalysis/iProFun_sim_20201228/sim", simulationRunLabels_char
)
# walk(.x = simAnalysisFolders, .f = dir.create)

# END Setup



######  Simulation Design  ####################################################
###  Import Data  ###
# Treated values
simDesign_ls <- readRDS(
  paste0(
    simResFolders[2], "/",
    dataFiles_char[
      str_detect(dataFiles_char, "indicatorMatrices")
    ]
  )
)

# Predictor Data Sets
dataFiles_char <-
  dataFiles_char[!str_detect(dataFiles_char, "indicatorMatrices")]


# Design points
designs_char <- unique(
  str_split_fixed(dataFiles_char, pattern = "_", n = 2)[, 2]
)

designGroups_ls <- 
  map(
    .x = designs_char,
    .f = ~{
      dataFiles_char[str_detect(dataFiles_char, pattern = .x)]
    }
  )

rm(dataFiles_char)

# Import the Three Data Sets with Strongest Signal
dataDF_ls <- map(
  .x = designGroups_ls[[19]],
  .f = ~{
    readRDS(paste0(simResFolders[2], "/", .x))
  }
)
names(dataDF_ls) <- c("CNV", "Prot", "RNAseq")

dataMatT_ls <- 
  dataDF_ls %>% 
  map(TransposeAssay) %>% 
  map(rename, Gene = Sample)



######  Analysis Internals  ###################################################
# We are using the defaults from their README
t0 <- Sys.time()
coad_iProFunPerm <- iProFun_permutate(
  ylist = dataMatT_ls[c("RNAseq", "Prot")],
  xlist = dataMatT_ls["CNV"],
  covariates = NULL,
  pi = rep(0.05, 2), # originally 3
  permutate_number = 10,
  grids = c(
    seq(0.75, 0.99, 0.01),
    seq(0.991, 0.999, 0.001),
    seq(0.9991, 0.9999, 0.0001)
  ),
  filter = 1,
  seed = 123
)
t1 <- Sys.time()

# 67.97 sec for one permutation; 8.124 min for 10 permutations; 70.84 min for 100
sigGenes_df <-
  coad_iProFunPerm$Gene_fdr[[1]] %>% 
  data.frame(stringsAsFactors = FALSE) %>% 
  rename(Gene = V1) %>% 
  rowwise() %>% 
  mutate(Signif = Y1 == "1" | Y2 == "1") %>% 
  ungroup() %>% 
  select(Gene, Signif)

iProFun_res <- 
  sigGenes_df %>% 
  filter(Signif) %>% 
  pull(Gene)
# The number of significant genes has not changed from 1 permutation, 10
#   permutations, or even 100 permutations. This is a problem.



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
  features = iProFun_res,
  universe = sigGenes_df$Gene
)


# Tabulate Results
iProFunRes_df <- 
  tibble(
    Pathway = cluster_PC$TERMS,
    FisherP = PW_p
  ) %>% 
  mutate(
    treated = Pathway %in% simDesign_ls$treatedPathways
  )

iProFunRes_df  %>% 
  filter(treated) %>% 
  arrange(FisherP)
# 40% power

# AUC?
as.numeric(
  pROC::roc(
    response = iProFunRes_df$treated,
    predictor = iProFunRes_df$FisherP,
    levels = as.factor(c(TRUE, FALSE)),
    # Case p-values should be lower than control p-values
    direction = "<"
  )$auc
)



######  Make a Function  ######################################################
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

RunIProFun <- function(simRunLabel, PC) {
  # browser()
  
  ###  1. Import Fixed Data  ###
  simResFolder <- paste0(
    "simResultsRaw/pathwayPCA_sim_20201201/sim", simRunLabel
  )
  # List all the data files for this run, and remove the design information that
  #   is stored in "indicatorMatricesXXX_ls.RDS".
  dataFiles_char <- list.files(simResFolder)
  
  # Treated values
  simDesign_ls <- readRDS(
    paste0(
      simResFolder, "/",
      dataFiles_char[
        str_detect(dataFiles_char, "indicatorMatrices")
      ]
    )
  )
  
  dataFiles_char <-
    dataFiles_char[!str_detect(dataFiles_char, "indicatorMatrices")]
  
  # Design points
  designs_char <- unique(
    str_split_fixed(dataFiles_char, pattern = "_", n = 2)[, 2]
  )
  
  designGroups_ls <- 
    map(
      .x = designs_char,
      .f = ~{
        dataFiles_char[str_detect(dataFiles_char, pattern = .x)]
      }
    )
  
  rm(dataFiles_char)
  
  
  ###  2. Loop Over Simulated Data  ###
  walk(
    .x = seq_along(designGroups_ls),
    .f = ~{
      # browser()
      
      ##  Data  ---
      dataDF_ls <- lapply(
        X = designGroups_ls[[.x]],
        FUN = function(fileName){ readRDS(paste0(simResFolder, "/", fileName)) }
      )
      names(dataDF_ls) <- c("CNV", "Prot", "RNAseq")
      dataMatT_ls <- 
        dataDF_ls %>% 
        map(TransposeAssay) %>% 
        map(rename, Gene = Sample)
      rm(dataDF_ls)
      
      
      ##  Analysis  ---
      t0 <- Sys.time()
      out <- iProFun_permutate(
        ylist = dataMatT_ls[c("RNAseq", "Prot")],
        xlist = dataMatT_ls["CNV"],
        covariates = NULL,
        pi = rep(0.05, 2), # originally 3
        permutate_number = 10,
        grids = c(
          seq(0.75, 0.99, 0.01),
          seq(0.991, 0.999, 0.001),
          seq(0.9991, 0.9999, 0.0001)
        ),
        filter = 1,
        seed = 123
      )
      t1 <- Sys.time()
      
      sigGenes_df <-
        out$Gene_fdr[[1]] %>% 
        data.frame(stringsAsFactors = FALSE) %>% 
        rename(Gene = V1) %>% 
        rowwise() %>% 
        mutate(Signif = Y1 == "1" | Y2 == "1") %>% 
        ungroup() %>% 
        select(Gene, Signif)
      
      
      ##  Wrangle Results  ---
      iProFun_res <- 
        sigGenes_df %>% 
        filter(Signif) %>% 
        pull(Gene)
      PW_p <- sapply(
        PC$pathways,
        testOverrep,
        features = iProFun_res,
        universe = sigGenes_df$Gene
      )
      # Tabulate Results
      iProFunRes_df <- 
        tibble(
          Pathway = PC$TERMS,
          FisherP = PW_p
        ) %>% 
        mutate(
          treated = Pathway %in% simDesign_ls$treatedPathways
        )
      auc_num <- as.numeric(
        pROC::roc(
          response = iProFunRes_df$treated,
          predictor = iProFunRes_df$FisherP,
          levels = as.factor(c(TRUE, FALSE)),
          # Case p-values should be lower than control p-values
          direction = "<"
        )$auc
      )
      
      
      ###  Save  ###
      results_ls <- list(
        Run = simResFolder,
        iProFun_results = out,
        enrichGenes = iProFun_res,
        pathwayPvals_df = iProFunRes_df,
        AUC = auc_num,
        computeTime = t1 - t0
      )
      
      cat("Saving iProFun results", designs_char[.x], "for run", simRunLabel, "at", format(Sys.time()), "\n")
      saveRDS(
        object = results_ls,
        file = paste0(
          "./simAnalysis/iProFun_sim_20201228/sim", simRunLabel, "/",
          "composite_", designs_char[.x]
        )
      )
      
    }
  )
  
}

# Test
RunIProFun("001", PC = cluster_PC)




######  Apply the Function  ###################################################
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

system.time(
  walk(
    .x = simulationRunLabels_char[2:100],
    RunIProFun,
    PC = cluster_PC
  )
)
# roughly 2 hrs and 26 minutes for one run. Based on this, we expect this to be
#   complete in roughly 10 days' time.
# 8.986668 days


