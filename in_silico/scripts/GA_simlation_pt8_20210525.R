# Calculate Single-Omics Results
# Gabriel Odom
# 2021-05-25

# In the last simulation step (scripts/GA_simulation_pt7_20201202.R), we
#   compiled the simulation results from 100 simulation replicates across 20
#   design points. Now, we would like to see (for ourselves), how well a single
#   -omic analysis would have worked.

library(tidyverse)
resultsPath <- "./simAnalysis/pathwayPCA_sim_20201202/"
simResults_ls <- readRDS(
  file = paste0(resultsPath, "sim_results_raw_20201202.RDS")
)

# END intro



######  Function to Calculate Single -Omics AUC and Test Size  ################
library(pROC)
SingleOmicsAUC <- function(run_df, platform = c("CNV", "RNAseq", "Prot")){
  # browser()
  
  platform <- match.arg(platform)
  platform <- paste0("rawp_", platform)
  
  AUC_num <- as.numeric(
    roc(
      response = run_df$treated,
      predictor = run_df[[platform]],
      levels = as.factor(c(TRUE, FALSE)),
      # Case p-values should be lower than control p-values
      direction = "<"
    )$auc
  )
  
  # This bit of witchcraft is from here:
  # https://shipt.tech/https-shipt-tech-advanced-programming-and-non-standard-evaluation-with-dplyr-e043f89deb3d
  signif_quo <- quo(!!sym(platform) < 0.05)
  testSize_num <- run_df %>%
    filter(!treated) %>%
    mutate(T1err = !!signif_quo) %>%
    pull(T1err) %>%
    mean()
  
  run_df %>% 
    select(Simulation, Partition, EffectSize) %>% 
    distinct() %>% 
    mutate(AUC = AUC_num, TSize = testSize_num)
  
}

# Test
SingleOmicsAUC(simResults_ls[[16]], platform = "Prot")



######  Apply, Wrangle, and Save the Results  #################################
AUC_df <- 
  map_dfr(
    .x = simResults_ls,
    .f = SingleOmicsAUC,
    platform = "Prot"
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
  # mutate(Method = "CNV Alone")
  # mutate(Method = "RNAseq Alone")
  mutate(Method = "Proteomics Alone")

saveRDS(
  aucClean_df,
  file = paste0(resultsPath, "SingleOmics_Prot_AUC_20210525.RDS")
)

# Plot these results in reports/simulation_v2_overview_20210330.Rmd

