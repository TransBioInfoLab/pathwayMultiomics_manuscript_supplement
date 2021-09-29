# Calculate MiniMax Power
# Gabriel Odom
# 2020-12-02

# We have 2 simulation replicates complete for pathwayPCA. I need to now write
#   a script that compiles power results from these data sets.
# NOTE: these results depend on the summary file created in this script when ran
#   on the FIU PC, because that is the computer that has all that data synced
#   locally. 

library(tidyverse)
resultsPath <- "./simAnalysis/pathwayPCA_sim_20201202/"

# END intro


######  Setup  ################################################################
simulationRunLabels_char <-
  seq_len(100) %>%
  as.character() %>%
  str_pad(width = 3, pad = "0")

simAnalysisFolders <- paste0(
  "simAnalysis/pathwayPCA_sim_20201202/sim", simulationRunLabels_char
)


# As of 2020-12-02, we only have the first 2 simulation replicates done.
# 60 complete on 2020-12-04; all by 2020-12-05
completeSims_idx <- 1:100

resultsFiles_ls <-
  map(
    .x = simAnalysisFolders[completeSims_idx],
    .f = ~{
      tibble(fileName = list.files(paste0(.x, "/")))
    }
  )
names(resultsFiles_ls) <- simulationRunLabels_char[completeSims_idx]

resultsFiles_df <-
  map_dfr(
    .x = resultsFiles_ls,
    .f = ~{
      .x %>%
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
  platforms_char <- c("CNV", "RNAseq", "Prot")
  resFiles_char <- paste0(
    platforms_char, "_",
    designRow_df[["Partition"]], "_",
    designRow_df[["EffectSize"]],
    ".RDS"
  )

  results_ls <-
    map(
      .x = paste0(resDir_path, resFiles_char),
      .f = readRDS
    ) %>%
    map("pVals_df") %>%
    map(., select, -FDR_BH)
  names(results_ls) <- platforms_char


  ###  Wrangle  ###
  results_ls %>%
    reduce(
      .f = left_join,
      by = c("terms", "treated")
    ) %>%
    rename(
      rawp_CNV    = rawp.x,
      rawp_RNAseq = rawp.y,
      rawp_Prot   = rawp
    ) %>%
    mutate(
      Simulation = designRow_df[["Simulation"]],
      Partition  = designRow_df[["Partition"]],
      EffectSize = designRow_df[["EffectSize"]]
    ) %>%
    select(
      Simulation, Partition, EffectSize, terms, treated, everything()
    )

}


# Test
View(
  ImportPowerResults(
    designRow_df = resultsFiles_df[1, ],
    res_dir = resultsPath
  )
)


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


# saveRDS(
#   simResults_ls,
#   file = paste0(resultsPath, "sim_results_raw_20201202.RDS")
# )
# 
# rm(list = ls())



######  Calculate MiniMax p-Value  ############################################
resultsPath <- "./simAnalysis/pathwayPCA_sim_20201202/"
simResults_ls <- readRDS(
  file = paste0(resultsPath, "sim_results_raw_20201202.RDS")
)

# betaParams_num <- c(alpha = 2, beta = 2)
# See "GA_simlation_pt6_20201202.R" for the ML estimation of alpha and beta
#   under H0
betaParams_num <- c(alpha = 1.847363, beta = 1.90439)

safeBetaMLE <- safely(Rfast::beta.mle)

MiniMaxTest <- function(df, alpha = 2, beta = 2, estimate = FALSE, ...){
  # browser()
  
  res1_df <- 
    df %>% 
    rowwise() %>% 
    mutate(
      MiniMaxTest = sort(
        c(rawp_CNV, rawp_RNAseq, rawp_Prot),
        decreasing = FALSE
      )[2]
    ) 
  
  # UPDATE 2021-05-14: include a way to estimate the parameters within each 
  #   simulation run.
  if(estimate) {
    # browser()
    
    nullTestStat_num <- 
      res1_df %>% 
      filter(!treated) %>% 
      pull(MiniMaxTest)
    
    beta_fit <- safeBetaMLE(nullTestStat_num, ...)
    
    if (!is.null(beta_fit$error)) {
      
      warning("MLE not convergent. Using supplied default values.")
      bestAlpha <- alpha
      bestBeta <- beta
      
    } else {
      bestAlpha <- beta_fit$result$param[["alpha"]]
      bestBeta <- beta_fit$result$param[["beta"]]
    }
    
  } else {
    bestAlpha <- alpha
    bestBeta <- beta
  }
  
  res1_df %>% 
    mutate(
      MiniMaxP = pbeta(q = MiniMaxTest, shape1 = bestAlpha, shape2 = bestBeta)
    ) %>% 
    ungroup()
  
}

# Test
test1 <- MiniMaxTest(
  df = simResults_ls[[1]], alpha = betaParams_num[1], beta = betaParams_num[2]
)
test2 <- MiniMaxTest(df = simResults_ls[[1]], alpha = 2, beta = 2)

View(MiniMaxTest(df = simResults_ls[[1]], estimate = TRUE))

# Apply (MLEs)
simMiniMax_ls <- 
  map(
    .x = simResults_ls,
    .f = MiniMaxTest,
    alpha = betaParams_num["alpha"],
    beta  = betaParams_num["beta"]
    # estimate = TRUE
  )



######  Calculate MiniMax Sensitivity and Specificity  ####################
###  AUC  ###
library(pROC)
MiniMaxAUC <- function(run_df){
  
  AUC_num <- as.numeric(
    roc(
      response = run_df$treated,
      predictor = run_df$MiniMaxP,
      # predictor = run_df$rawp_CNV,
      levels = as.factor(c(TRUE, FALSE)),
      # Case p-values should be lower than control p-values
      direction = "<"
    )$auc
  )
  
  testSize_num <- run_df %>%
    filter(!treated) %>%
    mutate(T1err = MiniMaxP < 0.05) %>%
    # mutate(T1err = rawp_CNV < 0.05) %>%
    pull(T1err) %>%
    mean()
  
  run_df %>% 
    select(Simulation, Partition, EffectSize) %>% 
    distinct() %>% 
    mutate(AUC = AUC_num, TSize = testSize_num)
  
}

# Test
MiniMaxAUC(simMiniMax_ls[[16]])

AUC_df <- 
  map_dfr(
    .x = simMiniMax_ls,
    .f = MiniMaxAUC
  )

# UPDATE 2021-05-21: 
# SINGLE-OMICS CODE MOVED TO scripts/GA_simulation_pt8_20210525.R


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
  # mutate(Method = "MiniMax (1.85,1.9)") # MiniMax_MLE_AUC_20210514.RDS
  # mutate(Method = "MiniMax (2,2)")      # MiniMax_2_2_AUC_20210513.RDS
  mutate(Method = "MiniMax (run est.)")   # MiniMax_MLEn45_AUC_20210514.RDS

# From Lily: we are calling this method MiniMax Test and MiniMax p-Value

saveRDS(
  aucClean_df,
  file = paste0(resultsPath, "MiniMax_MLEn45_AUC_20210514.RDS")
)



###  Power  ###
MiniMaxPower <- function(df){
  
  df %>% 
    group_by(Simulation, Partition, EffectSize) %>% 
    summarise(
      TP = sum(MiniMaxP <  0.05 &  treated),
      TN = sum(MiniMaxP >= 0.05 & !treated),
      FP = sum(MiniMaxP <  0.05 & !treated),
      FN = sum(MiniMaxP >= 0.05 &  treated),
      Sensitivity = TP / (FN + TP),
      Specificity = TN / (FP + TN),
      .groups = "drop"
    )
  
}

# Test
MiniMaxPower(df = simMiniMax_ls[[1]])
MiniMaxPower(df = simMiniMax_ls[[12]])
MiniMaxPower(df = simMiniMax_ls[[20]])

# Apply
power_df <- 
  map_dfr(
    .x = simMiniMax_ls,
    .f = MiniMaxPower
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

# From Lily: we are calling this method MiniMax Test and MiniMax p-Value

saveRDS(
  powerClean_df,
  file = paste0(resultsPath, "MiniMax_power_20201202.RDS")
)

rm(list = ls())



######  Plot AUC  ###########################################################
resultsPath <- "./simAnalysis/pathwayPCA_sim_20201202/"
miniMaxAUC_df <- readRDS(
  file = paste0(resultsPath, "MiniMax_AUC_20201202.RDS")
)

ggplot(data = miniMaxAUC_df) +
  theme_bw() +
  aes(x = EffectSize, y = AUC, group = EffectSize) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(
    title = "MiniMax Detection AUC by Treatment Proportion and Effect Size",
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
resultsPath <- "./simAnalysis/pathwayPCA_sim_20201202/"
miniMaxPowerClean_df <- readRDS(
  file = paste0(resultsPath, "MiniMax_power_20201202.RDS")
)

ggplot(data = miniMaxPowerClean_df) +
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
