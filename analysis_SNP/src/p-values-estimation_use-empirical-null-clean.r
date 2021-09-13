setwd("C:/Users/lxw391/TBL Dropbox/Lily Wang/pathwayPCA_multiOmics/LWtest/GWAS_pathways")

# use empirical null distribution to estimate pathway p-values 
library(readxl)
data <- read_excel("C:/Users/lxw391/TBL Dropbox/Lily Wang/pathwayPCA_multiOmics/LWtest/GWAS_pathways/RESULT/all_results_7-16-2021.xlsx")

############## empirical null - bacon approach
library(bacon)

bc <- bacon(
  teststatistics = NULL,
  effectsizes =  data$Estimate,
  standarderrors = data$StdErr,
  na.exclude = TRUE
)

data.with.inflation <- data.frame(
  data,
  Estimate.bacon = bacon::es(bc),
  StdErr.bacon = bacon::se(bc),
  pValue.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE
)

#### inflation factor - after bacon correction

bc <- bacon(
  teststatistics = NULL,
  effectsizes =  data.with.inflation$Estimate.bacon,
  standarderrors = data.with.inflation$StdErr.bacon,
  na.exclude = TRUE
)
# inflation factor
print("lambda.bacon")
print(inflation(bc))   #0.9402

# traditional lambda
data.with.inflation$zvalue.bacon <- data.with.inflation$Estimate.bacon / data.with.inflation$StdErr.bacon
data.with.inflation$chisq.bacon <- (data.with.inflation$zvalue.bacon) ^ 2
inflationFactor <- median(data.with.inflation$chisq.bacon,na.rm = TRUE) / qchisq(0.5, 1)
print("lambda")
print(inflationFactor)  #0.956

write.csv (data.with.inflation, "./RESULT/all_results_bacon_correction_7-21-2021.csv", row.names = FALSE)

