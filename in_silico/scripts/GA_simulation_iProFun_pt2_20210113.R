# Calculate iProFun:: Power
# Gabriel Odom
# 2021-01-13
# UPDATED: 2021-05-02

# We have all simulation replicates complete for iProFun. I need to now write
#   a script that compiles power results from these data sets.
# NOTE: these results depend on the summary file created in this script when ran
#   on the FIU PC, because that is the computer that has all that data synced
#   locally. 

library(tidyverse)
resultsPath <- "./simAnalysis/iProFun_sim_20201228/"

# END intro


######  Setup  ################################################################
simulationRunLabels_char <-
  seq_len(100) %>%
  as.character() %>%
  str_pad(width = 3, pad = "0")

simAnalysisFolders <- paste0(
  "simAnalysis/iProFun_sim_20201228/sim", simulationRunLabels_char
)

resultsFiles_ls <-
  map(
    .x = simAnalysisFolders,
    .f = ~{ list.files(paste0(.x, "/")) }
  )
names(resultsFiles_ls) <- simulationRunLabels_char

resultsFiles_df <-
  map_dfr(
    .x = resultsFiles_ls,
    .f = ~{
      tibble(fileName = .x) %>%
        mutate(fileName = str_remove(fileName, ".RDS")) %>%
        separate(
          fileName,
          into = c("Type", "Partition", "EffectSize"),
          sep = "_"
        ) %>%
        select(-Type) %>%
        distinct()
    },
    .id = "Simulation"
  )



######  Import Power Results  #################################################
ImportPowerResults <- function(designRow_df, res_dir){
  # browser()
  
  ###  Import  ###
  resDir_path <- paste0(
    res_dir, "sim", designRow_df[["Simulation"]], "/"
  )
  resFile_char <- paste0(
    "composite_",
    designRow_df[["Partition"]], "_",
    designRow_df[["EffectSize"]], ".RDS"
  )
  
  results_ls <- readRDS(paste0(resDir_path, resFile_char))
  
  
  ###  Wrangle  ###
  results_ls %>%
    pluck("pathwayPvals_df") %>%
    mutate(
      Simulation = designRow_df[["Simulation"]],
      Partition  = designRow_df[["Partition"]],
      EffectSize = designRow_df[["EffectSize"]]
    ) %>%
    select(
      Simulation, Partition, EffectSize, everything()
    )
  
}


# # Test
# View(
#   ImportPowerResults(
#     designRow_df = resultsFiles_df[40, ],
#     res_dir = resultsPath
#   )
# )


# Apply
simResults_ls <-
  map(
    .x = seq_len(nrow(resultsFiles_df)),
    .f = ~{
      ImportPowerResults(
        designRow_df = resultsFiles_df[.x, ],
        res_dir = resultsPath
      )
    }
  )


saveRDS(
  object = simResults_ls,
  file = paste0(resultsPath, "sim_results_raw_20210113.RDS")
)

rm(list = ls())




######  Calculate MiniMax Sensitivity and Specificity  ####################
resultsPath <- "./simAnalysis/iProFun_sim_20201228/"
simResults_ls <- readRDS(
  file = paste0(resultsPath, "sim_results_raw_20210113.RDS")
)


###  AUC  ###
library(pROC)
iProFunAUC <- function(run_df){
  
  AUC_num <- as.numeric(
    roc(
      response = run_df$treated,
      predictor = run_df$FisherP,
      levels = as.factor(c(TRUE, FALSE)),
      # Case p-values should be lower than control p-values
      direction = "<"
    )$auc
  )
  
  testSize_num <- run_df %>%
    filter(!treated) %>%
    mutate(T1err = FisherP < 0.05) %>%
    pull(T1err) %>%
    mean()
  
  run_df %>% 
    select(Simulation, Partition, EffectSize) %>% 
    distinct() %>% 
    mutate(AUC = AUC_num, TSize = testSize_num)
  
}

# Test
iProFunAUC(simResults_ls[[60]])

AUC_df <- 
  map_dfr(
    .x = simResults_ls,
    .f = iProFunAUC
  )

aucClean_df <- 
  AUC_df %>% 
  mutate(
    Partition  = str_remove(Partition, "partition"),
    EffectSize = str_remove(EffectSize, "delta")
  ) %>% 
  mutate(
    # UPDATE 2020-09-01: We changed from tenths to fifths
    TreatProp  = as.numeric(Partition) / 5,
    EffectSize = as.numeric(EffectSize)
  ) %>% 
  select(-Partition) %>% 
  select(Simulation, TreatProp, EffectSize, everything()) %>% 
  mutate(Method = "iProFun")

saveRDS(
  aucClean_df,
  file = paste0(resultsPath, "iProFun_AUC_20210502.RDS")
)



###  Power  ###
iProFunPower <- function(df){
  # browser()
  
  df %>% 
    # mutate(FisherFDR = p.adjust(FisherP, method = "fdr")) %>% 
    group_by(Simulation, Partition, EffectSize) %>% 
    summarise(
      TP = sum(FisherP <  0.05 &  treated),
      TN = sum(FisherP >= 0.05 & !treated),
      FP = sum(FisherP <  0.05 & !treated),
      FN = sum(FisherP >= 0.05 &  treated),
      Sensitivity = TP / (FN + TP),
      Specificity = TN / (FP + TN),
      .groups = "drop"
    )
  
}

# Test
iProFunPower(df = simResults_ls[[1]])
iProFunPower(df = simResults_ls[[12]])
iProFunPower(df = simResults_ls[[20]])
iProFunPower(df = simResults_ls[[32]])
iProFunPower(df = simResults_ls[[40]])

# Apply
power_df <- 
  map_dfr(
    .x = simResults_ls,
    .f = iProFunPower
  )

powerClean_df <- 
  power_df %>% 
  mutate(
    Partition  = str_remove(Partition, "partition"),
    EffectSize = str_remove(EffectSize, "delta")
  ) %>% 
  mutate(
    # UPDATE 2020-09-01: We changed from tenths to fifths
    TreatProp  = as.numeric(Partition) / 5,
    EffectSize = as.numeric(EffectSize)
  ) %>% 
  select(-Partition) %>% 
  select(Simulation, TreatProp, EffectSize, everything())

saveRDS(
  powerClean_df,
  file = paste0(resultsPath, "iProFun_power_20210116.RDS")
)

rm(list = ls())



######  Plot AUC  ###########################################################
resultsPath <- "./simAnalysis/iProFun_sim_20201228/"
iProFunAUC_df <- readRDS(
  file = paste0(resultsPath, "iProFun_AUC_20210116.RDS")
)

ggplot(data = iProFunAUC_df) +
  theme_bw() +
  aes(x = EffectSize, y = AUC, group = EffectSize) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(
    title = "iProFun Detection AUC by Treatment Proportion and Effect Size",
    subtitle = "Facets are the proportion of treated features in pathways.",
    x = "Increased Signal Strength (10% - 50%)",
    y = "Area Under ROC Curve (AUC)"
  ) +
  geom_violin() +
  geom_boxplot() +
  # geom_jitter(height = 0.01, alpha = 0.1, size = 2) +
  geom_hline(yintercept = 0.8, colour = "blue") +
  facet_wrap(~TreatProp)



######  Plot Power  ###########################################################
resultsPath <- "./simAnalysis/iProFun_sim_20201228/"
iProFunPower_df <- readRDS(
  file = paste0(resultsPath, "iProFun_power_20210116.RDS")
)

ggplot(data = iProFunPower_df) +
  theme_bw() +
  aes(x = EffectSize, y = Sensitivity, group = EffectSize) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "iProFun Detection Power by Treatment Proportion and Effect Size",
    subtitle = "Facets are the proportion of treated features in pathways.",
    x = "Increased Signal Strength (10% - 50%)",
    y = "Power"
  ) +
  geom_violin() +
  geom_boxplot() +
  # geom_jitter(height = 0.01, alpha = 0.1, size = 2) +
  geom_hline(yintercept = 0.8, colour = "blue") +
  facet_wrap(~TreatProp)

# Remember that power is meaningless without a well-controlled test size.
