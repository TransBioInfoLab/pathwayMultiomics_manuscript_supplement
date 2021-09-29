# Analyse Simulated Data with pathwayPCA
# Gabriel Odom
# 2020-08-26


library(tidyverse)
library(pathwayPCA)

cluster_PC <- read_gmt(file = "data/cluster_pathway_collection_20201117.gmt")


######  Directory Setup  ######################################################
simulationRunLabels_char <- 
  seq_len(100) %>% 
  as.character() %>% 
  str_pad(width = 3, pad = "0")

simResFolders <- paste0(
  "simResultsRaw/pathwayPCA_sim_20201201/sim", simulationRunLabels_char
)
dataFiles_char <- list.files(simResFolders[1])

simAnalysisFolders <- paste0(
  "simAnalysis/pathwayPCA_sim_20201202/sim", simulationRunLabels_char
)
# walk(.x = simAnalysisFolders, .f = dir.create)

# END Setup



######  Simulation Design  ####################################################
###  Import Data  ###
# Treated values
simDesign_ls <- readRDS(
  paste0(
    simResFolders[1], "/",
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

rm(dataFiles_char, designs_char)

# Import the Three Data Sets with Strongest Signal
dataDF_ls <- map(
  .x = designGroups_ls[[20]],
  .f = ~{
    readRDS(paste0(simResFolders[1], "/", .x))
  }
)
names(dataDF_ls) <- c("CNV", "Prot", "RNAseq")

# Extract response
# Which rows have treatment indicators?
treatedSamps_idx <- unique(simDesign_ls$indicatorMatrices$partition1$i)
resp_df <- tibble(
  Sample = simDesign_ls$indicatorDimNames[[1]],
  Treated = Sample %in% simDesign_ls$indicatorDimNames[[1]][treatedSamps_idx]
)



######  Inspect Treated Data  #################################################
# Can we see a treatment effect? If we can't see one, then perhaps this is why
#   the methods we have tried so far can't see one either.

# testPath_ls <- cluster_PC[[simDesign_ls$treatedPathways[1]]]
# 
# test_df <- dataDF_ls$Prot[, c("Sample", testPath_ls$IDs)]
# treatAve_num <- 
#   test_df %>%
#   filter(
#     Sample %in% pull(filter(resp_df, Treated), Sample)
#   ) %>%
#   select(-Sample) %>%
#   colMeans
# unTreatAve_num <- 
#   test_df %>%
#   filter(
#     Sample %in% pull(filter(resp_df, !Treated), Sample)
#   ) %>%
#   select(-Sample) %>%
#   colMeans
# 
# hist(treatAve_num - unTreatAve_num, breaks = 10)
# rm(testPath_ls, test_df, treatAve_num, unTreatAve_num)

# Recall that 80% of the genes in the selected pathway receive treatment, so it
#   makes sense that this is bimodal. There is a clear treatment effect, as
#   80% of these differences are > 0. the question becomes: can pathwayPCA + the
#   MiniMax statistic detect this difference? So far, sCCA, iProFun, and ADE4
#   have not been able to detect this difference.




######  Analysis Internals  ###################################################

###  1. Read Data and Create Omics*  ###
# List all the data files for this run, and remove the design information that
#   is stored in "indicatorMatricesXXX_ls.RDS".
omics_ls <- map(
  .x = dataDF_ls,
  .f = CreateOmics,
  pathwayCollection_ls = cluster_PC,
  response = resp_df,
  respType = "categ"
)

map(omics_ls, print)

###  2. AES-PCA and p-Values  ###
start_POSIX <- Sys.time()
res_aespcOut <- AESPCA_pVals(
  omics_ls[[1]],
  parallel = TRUE, numCores = 5L,
  adjustpValues = TRUE,
  adjustment = "BH"
)
end_POSIX <- Sys.time()
res_aespcOut$compTime <- end_POSIX - start_POSIX
# 1.557 min for 2 cores; 1.117 min for 4; 1.056 for 6; 1.120 for 8; 0.970 for 5


res_aespcOut$pVals_df <-
  res_aespcOut$pVals_df %>% 
  mutate(treated = terms %in% simDesign_ls$treatedPathways) %>% 
  select(terms, rawp, FDR_BH, treated) 

saveRDS(
  res_aespcOut, paste0(simAnalysisFolders[78], "/", allSimDataFiles_char[60])
)



######  Make a Function  ######################################################
AnalyzeSimulation <- function(simDataFolder, simAnalFolder, PC){
  # browser()
  
  ###  1. Import Fixed Data  ###
  # List all the data files for this run, and remove the design information that
  #   is stored in "indicatorMatricesXXX_ls.RDS".
  allSimDataFiles_char <- list.files(simDataFolder)
  
  # Treated values
  simDesign_ls <- readRDS(
    paste0(
      simDataFolder, "/",
      allSimDataFiles_char[
        str_detect(allSimDataFiles_char, "indicatorMatrices")
      ]
    )
  )
  # Which rows have treatment indicators?
  treatedSamps_idx <- unique(simDesign_ls$indicatorMatrices[[1]]$i)
  resp_df <- tibble(
    Sample = simDesign_ls$indicatorDimNames[[1]],
    Treated = Sample %in% simDesign_ls$indicatorDimNames[[1]][treatedSamps_idx]
  )
  
  # Predictor Data Sets
  allSimDataFiles_char <-
    allSimDataFiles_char[
      !str_detect(allSimDataFiles_char, "indicatorMatrices")
    ]
  
  
  ###  2. Loop Over Simulated Data  ###
  walk(
    .x = allSimDataFiles_char,
    .f = ~{
      # browser()
      
      # Data
      sim_df <- readRDS(paste0(simDataFolder, "/", .x))
      omics <- CreateOmics(
        assayData_df = sim_df,
        pathwayCollection_ls = PC,
        response = resp_df,
        respType = "categ"
      )
      
      
      # Analysis
      start_POSIX <- Sys.time()
      res_aespcOut <- AESPCA_pVals(
        omics,
        parallel = TRUE, numCores = 5L,
        adjustpValues = TRUE,
        adjustment = "BH"
      )
      end_POSIX <- Sys.time()
      
      res_aespcOut$compTime <- end_POSIX - start_POSIX
      res_aespcOut$pVals_df <-
        res_aespcOut$pVals_df %>% 
        mutate(treated = terms %in% simDesign_ls$treatedPathways) %>% 
        select(terms, rawp, FDR_BH, treated) 
      
      
      # Save
      saveRDS(
        res_aespcOut, paste0(simAnalFolder, "/", .x)
      )
      
    }
  )
  
}

# # Test
# AnalyzeSimulation(
#   simDataFolder = simResFolders[1],
#   simAnalFolder = simAnalysisFolders[1],
#   PC = cluster_PC
# )



######  Apply the Function  ###################################################
t0 <- Sys.time()
walk2(
  .x = simResFolders[61:100],
  .y = simAnalysisFolders[61:100],
  AnalyzeSimulation,
  PC = cluster_PC
)
t1 <- Sys.time()
# 1.474341 hrs for first 2; 12.48115 for 3:20; about 5.5 hours for 21:28;
#   22.08815 hrs for 29:60; 27.87096 hrs for 61:100


