# Analyse Simulated Data with mitch::
# Gabriel Odom
# 2020-12-14

# It's taken a few months, but we have finally found a worthy competitor. The
#   mitch:: package
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06856-9
# https://github.com/markziemann/mitch
#
# This package has great power, but no way to control the test size. Thus, Lily
#   suggests that we use the AUC to compare this method. I've gone back and
#   added AUC calculations for the MiniMax (in GA_simlation_pt7_20201202.R and
#   MiniMax_wrangle_power_20200821.R). The MiniMax performs strongly with this 
#   new metric as well.
# Now that we have a comparable method and an appropriate metric, this script
#   will attempt to analyse all the simulated data with the mitch:: package.


# BiocManager::install("mitch")
library(tidyverse)
library(mitch)

clusterPC_ls <- gmt_import("data/cluster_pathway_collection_20201117.gmt")



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
  "simAnalysis/mitch_sim_20201214/sim", simulationRunLabels_char
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
  .x = designGroups_ls[[20]],
  .f = ~{
    readRDS(paste0(simResFolders[2], "/", .x))
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



######  Analysis Internals  ###################################################

start_POSIX <- Sys.time()
###  1. Single-Gene p-Values  ###
# List all the data files for this run, and remove the design information that
#   is stored in "indicatorMatricesXXX_ls.RDS".
metricValues_ls <- map(
  .x = dataDF_ls,
  .f = ~{
    # browser()
    
    data_df <- inner_join(resp_df, .x, by = "Sample")
    metric_num <- vector(mode = "double", length = ncol(.x) - 1)
    
    for(k in seq_along(metric_num)) { 
      
      geneExpr_num <- data_df[, c(k + 2), drop = TRUE]
      modSummary <- summary(lm(geneExpr_num ~ data_df$Treated))
      # Estimate: col 1; t-stat: col 3, p-value: col 4
      metric_num[k] <- coefficients(modSummary)[2, 3]
      
    }
    
    data.frame(metric_num)
    
  }
)

metricValues_df <- bind_cols(metricValues_ls)
colnames(metricValues_df) <- names(metricValues_ls)
rownames(metricValues_df) <- colnames(dataDF_ls$CNV[, -1])


###  2. mitch Composite p-Values  ###
out_mitch <- mitch_calc(
  x = metricValues_df,
  genesets = clusterPC_ls,
  minsetsize = 5,
  resrows = 10,
  priority = "effect"
)

end_POSIX <- Sys.time()
out_mitch$compTime <- end_POSIX - start_POSIX


###  3. Wrangle Results  ###
out_mitch$enrichment_result <-
  out_mitch$enrichment_result %>% 
  mutate(Treated = set %in% simDesign_ls$treatedPathways) 

saveRDS(
  out_mitch, paste0(simAnalysisFolders[1], "/composite_tStat_", designs_char[20])
)



######  Make a Function  ######################################################
AnalyzeSimulation <- function(simDataFolder, simAnalFolder, PC, metric = 4L){
  # metric is which column of the linear model coefficients matrix we care about
  #  Estimate: col 1; t-stat: col 3, p-value: col 4
  
  # browser()
  
  ###  1. Import Fixed Data  ###
  # List all the data files for this run, and remove the design information that
  #   is stored in "indicatorMatricesXXX_ls.RDS".
  dataFiles_char <- list.files(simDataFolder)
  
  # Treated values
  simDesign_ls <- readRDS(
    paste0(
      simDataFolder, "/",
      dataFiles_char[
        str_detect(dataFiles_char, "indicatorMatrices")
      ]
    )
  )
  # Which rows have treatment indicators?
  treatedSamps_idx <- unique(simDesign_ls$indicatorMatrices[[1]]$i)
  resp_df <- tibble(
    Sample = simDesign_ls$indicatorDimNames[[1]],
    Treated = Sample %in% simDesign_ls$indicatorDimNames[[1]][treatedSamps_idx]
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
        FUN = function(x){
          readRDS(paste0(simDataFolder, "/", x))
        }
      )
      names(dataDF_ls) <- c("CNV", "Prot", "RNAseq")
      
      
      ##  Analysis  ---
      start_POSIX <- Sys.time()
      
      # Single-Gene Testing
      metricValues_ls <- lapply(
        X = dataDF_ls,
        FUN = function(x){
          
          data_df <- inner_join(resp_df, x, by = "Sample")
          metric_num <- vector(mode = "double", length = ncol(x) - 1) + 1
          
          for(k in seq_along(metric_num)) { 
            
            geneExpr_num <- data_df[, c(k + 2), drop = TRUE]
            modSummary <- summary(lm(geneExpr_num ~ data_df$Treated))
            metric_num[k] <- coefficients(modSummary)[2, metric]
            
          }
          
          data.frame(metric_num)
          
        }
      )
      
      metricValues_df <- suppressMessages(bind_cols(metricValues_ls))
      colnames(metricValues_df) <- names(metricValues_ls)
      rownames(metricValues_df) <- colnames(dataDF_ls$CNV[, -1])
      
      # mitch
      out_mitch <- mitch_calc(
        x = metricValues_df,
        genesets = PC,
        minsetsize = 5,
        resrows = 10,
        priority = "effect"
      )
      
      end_POSIX <- Sys.time()
      out_mitch$compTime <- end_POSIX - start_POSIX
      
      
      ##  Wrangle and Save  ---
      out_mitch$enrichment_result <-
        out_mitch$enrichment_result %>% 
        mutate(Treated = set %in% simDesign_ls$treatedPathways) 
      
      cat("Saving mitch results", designs_char[.x], "at", format(Sys.time()), "\n")
      saveRDS(
        out_mitch,
        paste0(
          simAnalFolder, "/composite_metric", metric, "_", designs_char[.x]
        )
      )
      
    }
  )
  
}

# Test
AnalyzeSimulation(
  simDataFolder = simResFolders[2],
  simAnalFolder = simAnalysisFolders[2],
  PC = clusterPC_ls,
  metric = 3
)



######  Apply the Function  ###################################################
t0 <- Sys.time()
walk2(
  .x = simResFolders[11:100],
  .y = simAnalysisFolders[11:100],
  AnalyzeSimulation,
  PC = clusterPC_ls,
  metric = 3
)
t1 <- Sys.time()
# 3.850679 hrs for all 100; 12.67 min for first 5; 1.74148 hrs for 6:50.
# 28.2085 min for first 10, 3.474936 hrs for remaining 90
