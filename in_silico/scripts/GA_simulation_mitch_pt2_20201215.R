# Calculate mitch:: Power
# Gabriel Odom
# 2020-12-15
# UPDATE 2021-05-02

# We have all simulation replicates complete for mitch. I need to now write
#   a script that compiles power results from these data sets.
# NOTE: these results depend on the summary file created in this script when ran
#   on the FIU PC, because that is the computer that has all that data synced
#   locally. 
# DEFN: "metric 3" is the t-stat; "metric 4" or "" is the p-value

library(tidyverse)
resultsPath <- "./simAnalysis/mitch_sim_20201214/"

# END intro


######  Setup  ################################################################
simulationRunLabels_char <-
  seq_len(100) %>%
  as.character() %>%
  str_pad(width = 3, pad = "0")

simAnalysisFolders <- paste0(
  "simAnalysis/mitch_sim_20201214/sim", simulationRunLabels_char
)

resultsFiles_ls <-
  map(
    .x = simAnalysisFolders,
    .f = ~{
      tibble(fileName = list.files(paste0(.x, "/"))) %>% 
        mutate(tFiles = str_detect(fileName, pattern = "metric3")) %>% 
        filter(tFiles) %>% 
        select(-tFiles)
    }
  )
names(resultsFiles_ls) <- simulationRunLabels_char

resultsFiles_df <-
  map_dfr(
    .x = resultsFiles_ls,
    .f = ~{
      .x %>%
        mutate(fileName = str_remove(fileName, ".RDS")) %>%
        separate(
          fileName,
          into = c("Type", "Metric", "Partition", "EffectSize"),
          # into = c("Type", "Partition", "EffectSize"),
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
    designRow_df[["Metric"]], "_",
    designRow_df[["Partition"]], "_",
    designRow_df[["EffectSize"]], ".RDS"
  )
  
  results_ls <- readRDS(paste0(resDir_path, resFile_char))
  
  
  ###  Wrangle  ###
  results_ls %>%
    pluck("enrichment_result") %>%
    mutate(
      Simulation = designRow_df[["Simulation"]],
      Partition  = designRow_df[["Partition"]],
      EffectSize = designRow_df[["EffectSize"]]
    ) %>%
    select(
      Simulation, Partition, EffectSize, set, setSize, Treated, everything()
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
  file = paste0(resultsPath, "sim_results_metric3_raw_20201216.RDS")
)

rm(list = ls())




######  Calculate MiniMax Sensitivity and Specificity  ####################
resultsPath <- "./simAnalysis/mitch_sim_20201214/"
simResults_ls <- readRDS(
  # # t-Stat
  # file = paste0(resultsPath, "sim_results_metric3_raw_20201216.RDS")
  # # p-Value
  file = paste0(resultsPath, "sim_results_raw_20201215.RDS")
)


###  AUC  ###
library(pROC)
mitchAUC <- function(run_df){
  
  AUC_num <- as.numeric(
    roc(
      response = run_df$Treated,
      predictor = run_df$p.adjustMANOVA,
      levels = as.factor(c(TRUE, FALSE)),
      # Case p-values should be lower than control p-values
      direction = "<"
    )$auc
  )
  
  testSize_num <- run_df %>%
    filter(!Treated) %>%
    # The AUC is very high, but the test size is atrocious. This tells me that
    #   the 0.05 threshold is too high. I've turned it all the way down to 
    #   10^-12, and it's still a bit inflated (median = 0.15). The test size is
    #   under control at 10^-20 for the t-Stat metric, while the proper cutoff
    #   for the p-value metric is 10^-13.
    mutate(T1err = p.adjustMANOVA < (10 ^ -13)) %>%
    pull(T1err) %>%
    mean()
  
  run_df %>% 
    select(Simulation, Partition, EffectSize) %>% 
    distinct() %>% 
    mutate(AUC = AUC_num, TSize = testSize_num)
  
}

# Test
mitchAUC(simResults_ls[[40]])

AUC_df <- 
  map_dfr(
    .x = simResults_ls,
    .f = mitchAUC
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
  # mutate(Method = "mitch_tStat")
  mutate(Method = "mitch_pValue")

saveRDS(
  aucClean_df,
  # # t-Stat
  # file = paste0(resultsPath, "mitch_metric3_AUC_20210513.RDS")
  # # p-value
  file = paste0(resultsPath, "mitch_metric4_AUC_20210513.RDS")
)



###  Power  ###
mitchPower <- function(df){
  
  df %>% 
    group_by(Simulation, Partition, EffectSize) %>% 
    summarise(
      TP = sum(p.adjustMANOVA <  0.05 &  Treated),
      TN = sum(p.adjustMANOVA >= 0.05 & !Treated),
      FP = sum(p.adjustMANOVA <  0.05 & !Treated),
      FN = sum(p.adjustMANOVA >= 0.05 &  Treated),
      Sensitivity = TP / (FN + TP),
      Specificity = TN / (FP + TN),
      .groups = "drop"
    )
  
}

# Test
mitchPower(df = simResults_ls[[1]])
mitchPower(df = simResults_ls[[12]])
mitchPower(df = simResults_ls[[20]])
mitchPower(df = simResults_ls[[32]])
mitchPower(df = simResults_ls[[40]])

# Apply
power_df <- 
  map_dfr(
    .x = simResults_ls,
    .f = mitchPower
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
  file = paste0(resultsPath, "mitch_metric3_power_20201216.RDS")
)

rm(list = ls())



######  Plot AUC  ###########################################################
resultsPath <- "./simAnalysis/mitch_sim_20201214/"
mitchAUC_df <- readRDS(
  file = paste0(resultsPath, "mitch_metric3_AUC_20201216.RDS")
)

ggplot(data = mitchAUC_df) +
  theme_bw() +
  aes(x = EffectSize, y = AUC, group = EffectSize) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(
    title = "mitch Detection AUC by Treatment Proportion and Effect Size",
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
resultsPath <- "./simAnalysis/mitch_sim_20201214/"
mitchPower_df <- readRDS(
  file = paste0(resultsPath, "mitch_power_20201215.RDS")
)

ggplot(data = mitchPower_df) +
  theme_bw() +
  aes(x = EffectSize, y = Sensitivity, group = EffectSize) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "MiniMax Detection Power by Treatment Proportion and Effect Size",
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