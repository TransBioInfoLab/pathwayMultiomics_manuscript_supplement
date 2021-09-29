# MiniMax Statistic Distribution and Test Size
# Gabriel Odom
# 2020-12-02



######  Introduction  #########################################################
# We performed a simulation on the test size of the MiniMax statistic. The
#   script is GA_simlation_pt5_20201202.R and the results are saved in
#   simResultsRaw/pathwayPCA_sim_20201201/miniMax_test_size_20201201.RDS.
#   These results have the MiniMax statistic values under no biological signal
#   (based on randomly-assigned clinical outcome). Using these values, in this
#   script we can: 
#   1. estimate the distribution of the MiniMax statistic under H0; specifically
#      we already know from the literature that the MiniMax ~ Beta, but now we
#      find the MLEs for alpha and beta.
#   2. measure the test size of the MiniMax statistic

results_ls <- readRDS(
  "simResultsRaw/pathwayPCA_sim_20201201/miniMax_test_size_20201201.RDS"
)



######  Estimating Alpha and Beta  ############################################
library(tidyverse)
library(Rfast)

results_df <- bind_rows(results_ls)
MiniMaxClean_num <- results_df$MiniMax
# The Rfast::beta.mle() function can't handle p-values near the boundary of [0,1]
MiniMaxClean_num[MiniMaxClean_num == 0] <- 1e-08
MiniMaxClean_num[MiniMaxClean_num == 1] <- 1 - 1e-08
bestParams <- beta.mle(MiniMaxClean_num)$param
#    alpha     beta 
# 1.906676 1.860925 



######  Estimate the Test Size  ###############################################
map_dbl(
  .x = results_ls,
  .f = ~{
    
    pVals <- pbeta(
      q = .x$MiniMax,
      shape1 = bestParams[1],
      shape2 = bestParams[2]
    )
    mean(pVals < 0.05)
    
  }
) %>% 
  summary()
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0200  0.0300  0.0498  0.0600  0.2600 
