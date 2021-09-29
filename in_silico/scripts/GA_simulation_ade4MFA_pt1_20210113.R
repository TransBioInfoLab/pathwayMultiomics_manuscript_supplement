# Analyse Simulated Data with MFA from ade4::
# Gabriel Odom
# 2021-01-13

# See the original script in "test_ade4_20201223.R"

# install.packages("ade4")
# BiocManager::install("made4")
library(tidyverse)
library(pathwayPCA)
library(ade4)

cluster_PC <- read_gmt(file = "data/cluster_pathway_collection_20201117.gmt")
allGenes_char <- 
  cluster_PC$pathways %>% 
  unlist() %>% 
  unique() %>% 
  sort()
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
  "simAnalysis/ade4_sim_20210113/sim", simulationRunLabels_char
)
# walk(.x = simAnalysisFolders, .f = dir.create)

# END Setup



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

RunMFA <- function(simRunLabel, PC, matchedGenes_char) {
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
        map(select, Sample, one_of(matchedGenes_char)) %>% 
        # ade4 requires "rivers" as rows
        map(TransposeAssay) %>% 
        map(tibble::column_to_rownames, var = "Sample")
      rm(dataDF_ls)
      
      
      ##  Analysis  ---
      t0 <- Sys.time()
      synth_ktab <- ktab.list.df(dataMatT_ls)
      out <- mfa(synth_ktab, option = "lambda1", scannf = FALSE, nf = 3)
      t1 <- Sys.time()
      
      geneP <- pnorm(
        out$l1[, 1],
        sd = 1 / sqrt(length(matchedGenes_char)),
        lower.tail = FALSE
      )
      
      
      ##  Wrangle Results  ---
      sigGenes_df <- tibble(
        Gene   = rownames(out$l1),
        Score  = out$l1[, 1],
        pVal   = geneP,
        Signif = geneP < 0.05
      )
      
      mfa_res <- 
        sigGenes_df %>% 
        filter(Signif) %>% 
        pull(Gene)
      PW_p <- sapply(
        PC$pathways,
        testOverrep,
        features = mfa_res,
        universe = sigGenes_df$Gene
      )
      # Tabulate Results
      mfaRes_df <- 
        tibble(
          Pathway = PC$TERMS,
          FisherP = PW_p
        ) %>% 
        mutate(
          treated = Pathway %in% simDesign_ls$treatedPathways
        )
      auc_num <- as.numeric(
        pROC::roc(
          response = mfaRes_df$treated,
          predictor = mfaRes_df$FisherP,
          levels = as.factor(c(TRUE, FALSE)),
          # Case p-values should be lower than control p-values
          direction = "<"
        )$auc
      )
      
      
      ###  Save  ###
      results_ls <- list(
        Run = simResFolder,
        mfa_results = out,
        enrichGenes = mfa_res,
        pathwayPvals_df = mfaRes_df,
        AUC = auc_num,
        computeTime = t1 - t0
      )
      
      cat("Saving MFA results", designs_char[.x], "for run", simRunLabel, "at", format(Sys.time()), "\n")
      saveRDS(
        object = results_ls,
        file = paste0(
          "./simAnalysis/ade4_sim_20210113/sim", simRunLabel, "/",
          "composite_", designs_char[.x]
        )
      )
      
    }
  )
  
}

# Test
RunMFA("001", PC = cluster_PC, matchedGenes_char = allGenes_char)




######  Apply the Function  ###################################################
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

system.time(
  walk(
    .x = simulationRunLabels_char[6:100],
    RunMFA,
    PC = cluster_PC,
    matchedGenes_char = allGenes_char
  )
)
# 15.146 min for first 5. We expect 4.8 hrs for the next 95; actual = 4.55 hrs.


