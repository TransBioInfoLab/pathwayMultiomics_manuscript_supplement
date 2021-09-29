# MiniMax Statistic Distribution and Test Size
# Gabriel Odom
# 2020-12-02



######  Introduction  #########################################################
# We performed a simulation on the test size of the MiniMax statistic. The
#   script is GA_simlation_pt5_20201202.R and the results are saved in
#   simResultsRaw/pathwayPCA_sim_20201201/globalTest_test_size_20201201.RDS.
#   These results have the MiniMax statistic values under no biological signal
#   (based on randomly-assigned clinical outcome). Using these values, in this
#   script we can: 
#   1. estimate the distribution of the MiniMax statistic under H0; specifically
#      we already know from the literature that the MiniMax ~ Beta, but now we
#      find the MLEs for alpha and beta.
#   2. measure the test size of the MiniMax statistic

results_ls <- readRDS(
  "simResultsRaw/pathwayPCA_sim_20201201/globalTest_test_size_20201201.RDS"
)



######  Estimating Alpha and Beta  ############################################
###  ML Estimates of alpha and beta Directly  ###
library(tidyverse)
library(Rfast)


betaMLE_ls <- map(
  .x = results_ls,
  .f = ~{
    
    beta_fit <- safely(beta.mle)(.x$MiniMax)
    
    if (!is.null(beta_fit$error)) {
      NULL
    } else {
      bestParams <- beta_fit$result
    }
    
  }
)

# What's our problem child?
map_lgl(betaMLE_ls, is.null) %>% which()
View(results_ls[[28]])
# Ah, I'll wager that the likelihood has an issue with boundary values. Let's
# test that:
beta.mle(runif(100))
beta.mle(c(0, runif(99)))
beta.mle(c(1, runif(99)))
# Yep, we can't have p-values on the boundary of [0, 1].

betaMLE_ls <- map(
  .x = results_ls,
  .f = ~{
    
    cleanX <- .x$MiniMax
    # Rfast's algorithm can't handle the boundaries, so use 10 * tolerance value
    #   from beta.mle()
    cleanX[cleanX == 0] <- 1e-08
    cleanX[cleanX == 1] <- 1 - 1e-08
    
    beta.mle(cleanX)$param
    
  }
)

map2_dbl(
  .x = results_ls,
  .y = betaMLE_ls,
  .f = ~{
    
    pVals <- pbeta(q = .x$MiniMax, shape1 = .y[1], shape2 = .y[2])
    mean(pVals < 0.05)
    
  }
) %>% 
  summary()
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0400  0.0400  0.0494  0.0600  0.1000 

# What do these parameters look like?
map_dfr(
  .x = betaMLE_ls,
  .f = ~{
    tibble(alpha = .x[1], beta = .x[2])
  }
) %>% 
  ggplot() +
  aes(x = alpha, y = beta) +
  geom_point()
# Not as tightly correlated as they were for the other simulations, but still
#   positive. Interesting.

map_dfr(
  .x = betaMLE_ls,
  .f = ~{
    tibble(alpha = .x[1], beta = .x[2])
  }
) %>% 
  colMeans()
# These averages are larger than 2, rather than slightly less than 2


###  Predicting the MLEs  ###
# We found (in "MiniMax_explore_distribution_w_correlation_20201105.R") that
#   there is a strong relationship between the correlation of the p-values and
#   the MLEs that best fit the Beta parameters alpha and beta. That relationship
#   was:
#      alpha = beta ~= (-7 / 8) * sqrt(1 - detCorr) + 2,
#   where detCorr is the determinant of the correlation matrix of the p-values.
map2_dfr(
  .x = results_ls,
  .y = betaMLE_ls,
  .f = ~{
    
    detCorr <- det(cor(.x[c("CNV_p", "RNA_p", "Prot_p")]))
    
    tibble(
      alpha = .y[1],
      beta = .y[2],
      predParam = (-7 / 8) * sqrt(1 - detCorr) + 2
    )
    
  }
)

# Well, this doesn't appear to work. In our simulations, we had 10,000 p-value
#   triplets which we used to estimate the MiniMax statistic. I believe that
#   this variance we currently see is due to the fact that we are estimating the
#   MiniMax statistic's distribution with only 50 p-value triplets.

# What if we concatenate all p-values and estimate the parameters that way?
results_df <- bind_rows(results_ls)
detCorr <- det(cor(results_df[c("CNV_p", "RNA_p", "Prot_p")]))

MiniMaxClean_num <- results_df$MiniMax
MiniMaxClean_num[MiniMaxClean_num == 0] <- 1e-08
MiniMaxClean_num[MiniMaxClean_num == 1] <- 1 - 1e-08
bestParams <- beta.mle(MiniMaxClean_num)$param
#    alpha     beta 
# 1.906676 1.860925 
(-7 / 8) * sqrt(1 - detCorr) + 2
# 1.887025
(1.860925 + 1.906676) / 2
# 1.8838

# And *now* it works. I guess we REALLY need the Law of Large Numbers for this
#   distribution to work.



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
# And we are good.
