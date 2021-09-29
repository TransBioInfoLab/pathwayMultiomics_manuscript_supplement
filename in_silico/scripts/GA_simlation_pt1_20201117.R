# Attempting Multi-Omics Simulation Data from TCGA Cases
# Gabriel Odom
# 2020-11-17


# UPDATE: this script is based on scripts/simulation_attempt_20200701.R. After
#   completing the simulation as we originally designed it, we realised that
#   none of the competing methods had comparable power (the three platforms as
#   is have only 5 samples in common). We will keep the results of this first
#   simulation, but Steven suggested that we also us the legacy format gene
#   expression data and repeat the simulation. Apparently, the "GenomeAnalyzer"
#   platform hs more overlap than the HiSeq platform (but it's a much older
#   technique, so the data has lower quality). Here's our updated plan:
#   1. Replace HiSeq RNAseq data with GA RNAseq data; this platform should have
#      nearly full overlap with the samples from the other platforms.
#   2. Repeat the MiniMax simulation study (technically, we should only have to
#      add a result for GA RNAseq because the MiniMax technique analyses the
#      platforms independently).
#   3. Once we have MiniMax results for a set of data with matched samples,
#      compare the MiniMax technique to some other methods (sCCA, NNMF, iProFun,
#      MoGSA, ade4, MSFA, etc.--I'll have to try these out). Obviously, for the
#      first simulation study (using Hi Seq RNAseq), the MiniMax statistic
#      offers vastly superior statistical power, because there are only 5 
#      samples available to the other methods. For the GenomeAnalyzer RNAseq,
#      there are closer to 220 matched samples, so we hope that the MiniMax will
#      not be too inferior.



######  Overview  #############################################################
# We're taking TCGA cases (from colon cancer) to simulate treated and untreated
#   samples. To set up the simulation study, we will
#   1. Cluster genes into 200 pathways with no overlap, p_1, ..., p_200. Pick
#      pathways with at least 5 genes in them (20% of 5 is 1 gene).
#   2. Identify cases (not necessarily matched) with CNV, protein, and RNAseq
#      values, n_CNV, n_prot, n_gene. Using cases will ensure that the
#      correlation structures will be the same for all data. Our technique
#      doesn't necessarily pick out second-moment differences.
#   3. Standardise all features to have mean 0 and standard deviation 1.
#   4. Ensure we have complete survival information for each subject (if we 
#      potentially would like to analyse survival response in addition to binary
#      classification).

# Given this pathway collection and these three data sets, we will simulate a
#   "subtype" of this cancer by treating half of the samples at random (ensuring
#   to treat those samples which are matched as matched. In order to do this, 
#   we should find the samples that are matched across pairs of the three data
#   sets and treat half of the samples in the intersection, then treat samples
#   outside these intersections until half of the samples from each data set
#   have been treated. This will require a table of treated and untreated sample
#   IDs for the three platforms for each simulation run.) In order to simulate
#   these data sets, we should use the simDataRaw/ directory. We can create
#   subdirectories by design point and sub-subdirectories by run ID. This will
#   have similar design to the "IamComparison" project of Pucher et al. (2019).

# Concerning the design points, we would like control over the following:
#   1. How many pathways of the 200 will be disregulated for the cancer subtype?
#      We will default to 20.
#      UPDATE 20201201: we have 50 pathways (because of the smaller set of
#      shared genes), and we will treat 5 of them.
#   2. How many of the pathways disregulated for the cancer will also be 
#      disregluated for each subject? (Not all pathways will be disregulated for
#      all subjects in real cancer.) At the start, this will default to 20, to
#      match the number of disregulated pathways; i.e., all pathways will be
#      disregulated for all subjects.
#      UPDATE 20201201: we've decided to hold off on this modification 
#   3. What proportion of genes in each disregulated pathway will be treated? We
#      will have this range over 20%, 40%, ..., 80%. 
#   4. What will the treatment effect be? Our first simulation with HiSeq RNAseq
#      had delta = +0.1x, +0.2x, ..., +0.5x of the standard deviation.


library(tidyverse)
library(pathwayPCA)



######  Colon Cancer Cases  ###################################################
# Grab data from the TCGA COADREAD page on:
# http://www.linkedomics.org/data_download/TCGA-COADREAD/

# RNAseq (HiSeq)
coadRNAseq_df <- read_delim(
  "data/Human__TCGA_COADREAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadRNAseq_df <- TransposeAssay(coadRNAseq_df)


# RNAseq (GenomeAnalyzer)
coadRNAseqGA_df <- read_delim(
  "data/Human__TCGA_COADREAD__UNC__RNAseq__GA_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadRNAseqGA_df <- TransposeAssay(coadRNAseqGA_df)


# Copy Number Variation (Log-Ratio)
coadCNV_df <- read_delim(
  "data/Human__TCGA_COADREAD__BI__SCNA__SNP_6.0__01_28_2016__BI__Gene__Firehose_GISTIC2.cct.gz", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadCNV_df <- TransposeAssay(coadCNV_df)


# Proteomics
coadProt_df <- read_delim(
  "data/Human__TCGA_COADREAD__VU__Proteome__Velos__01_28_2016__VU__Gene__CDAP_UnsharedPrecursorArea_r2.cct", 
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadProt_df <- TransposeAssay(coadProt_df)


# Clinical
coadClinical_df <- read_delim(
  "data/Human__TCGA_COADREAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
coadClinical_df <-
  TransposeAssay(coadClinical_df) %>% 
  select(Sample, overall_survival, status) %>% 
  drop_na()


###  Check Sample Overlaps  ###
# HiSeq RNAseq: 377; 367 after removing missing responses
length(intersect(coadRNAseq_df$Sample, coadClinical_df$Sample))
setdiff(coadRNAseq_df$Sample, coadClinical_df$Sample)

# GA RNAseq: 222; 205 after removing missing responses
length(intersect(coadRNAseqGA_df$Sample, coadClinical_df$Sample))
setdiff(coadRNAseqGA_df$Sample, coadClinical_df$Sample)

# CNV: 614; 582 after removing missing responses
length(intersect(coadCNV_df$Sample, coadClinical_df$Sample))
setdiff(coadCNV_df$Sample, coadClinical_df$Sample)

# Proteomics: 90; 81 after removing missing responses
length(intersect(coadProt_df$Sample, coadClinical_df$Sample))
setdiff(coadProt_df$Sample, coadClinical_df$Sample)
# Because we are going to simulate response anyway, we can use all the samples,
#   but we can use the matched samples to estimate the distribution of survival
#   times for the missing samples (this would be the best of both worlds). We
#   then state that the X genes chosen to be treated would be TSGs for our 
#   simulated cancer subtype, resulting in increased survival time (then we 
#   multiply survival time by 1.5 or 2 or something).


###  Check Gene Overlaps  ###
commonGenes_char <- 
  intersect(
    intersect(
      colnames(coadProt_df[, -1]),
      colnames(coadRNAseqGA_df[, -1])
    ),
    colnames(coadCNV_df[, -1])
  )
# There are 1710 genes shared by all three data sets
saveRDS(commonGenes_char, file = "data/common_genes_20201117.RDS")

length(
  unique(
    union_all(
      colnames(coadCNV_df[, -1]),
      colnames(coadRNAseqGA_df[, -1]),
      colnames(coadProt_df[, -1])
    )
  )
)

length(
  intersect(
    colnames(coadRNAseq_df[, -1]),
    colnames(coadRNAseqGA_df[, -1])
  )
)
# 5500 / 6149 genes shared between GA and HiSeq (large concordance, ~89%)

length(
  intersect(
    colnames(coadProt_df[, -1]),
    colnames(coadRNAseqGA_df[, -1])
  )
)
# GA: 1715 / 5538 ~= 31% concordance
# HiSeq: 5275 / 5538 ~ 95% concordance

length(
  intersect(
    colnames(coadCNV_df[, -1]),
    colnames(coadRNAseqGA_df[, -1])
  )
)
# 5212 / 6149 ~= 85%
# rm(coadClinical_df, coadCNV_df, coadProt_df, coadRNAseq_df)



######  Pathway Collection  ###################################################
# SEE scripts/GA_cluster_genes_into_pathways_20201117.R
cluster_PC <- read_gmt(file = "data/cluster_pathway_collection_20201117.gmt")



######  Sample Overlaps  ######################################################
allSamples_df <- tibble(
  Sample   = coadClinical_df$Sample,
  inRNAseq = Sample %in% coadRNAseqGA_df$Sample,
  inCNV    = Sample %in% coadCNV_df$Sample,
  inProt   = Sample %in% coadProt_df$Sample
)
allSamples_df %>% 
  select(-Sample) %>% 
  colSums()
# Checks out


###  Venn Diagram of Samples  ###
library(VennDiagram)
# See <https://www.r-graph-gallery.com/14-venn-diagramm.html> for help
library(RColorBrewer)
# display.brewer.all(colorblindFriendly = TRUE)
myCol3 <- brewer.pal(3, "Dark2")
myCol4 <- brewer.pal(4, "Dark2")

venn.diagram(
  # Base Plot
  x = list(
    coadCNV_df$Sample, coadRNAseqGA_df$Sample, coadProt_df$Sample
  ),
  category.names = c("CNV" , "RNAseq_GA" , "Prot"),
  filename = "results/tcga_coadread_samples_venn_20201117.png",
  output = TRUE,
  
  # Output features
  imagetype = "png" ,
  height = 1200, 
  width = 1600, 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  fill = myCol3,
  
  # Set names
  cat.default.pos = "outer",
  cat.pos = c(-25, 27, 180), # RNAseq is the first here
  cat.dist = c(0.055, 0.085, 0.085) # RNAseq is the second here
)

venn.diagram(
  # Base Plot
  x = list(
    coadCNV_df$Sample,
    coadRNAseqGA_df$Sample,
    coadProt_df$Sample,
    coadClinical_df$Sample
  ),
  category.names = c("CNV" , "RNAseq_GA" , "Prot", "Clinical"),
  filename = "results/tcga_coadread_samples_wClinical_venn_20201117.png",
  output = TRUE,
  
  # Output features
  imagetype = "png" ,
  height = 1200, 
  width = 1600, 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  fill = myCol4
)



######  Create Design Template  ###############################################
# Let's first compare to see if there is complete overlap between the data we've
#   already simulated and the new GA RNAseq platform (if there is, then we only
#   need to simulate new GA data). The simulation response map was created in
#   scripts/simulation_attempt_20200720.R.
simDesignTest_df <- read_csv(
  file = "simDataRaw/treated_samples_map_20200720.csv"
)

setdiff(coadRNAseqGA_df$Sample, unique(simDesignTest_df$Sample))
# All 222 samples from GA RNAseq have an assigned synthetic outcome
setdiff(unique(simDesignTest_df$Sample), coadRNAseqGA_df$Sample)
# 409 subjects are not included

# Based on this information, we should be able to keep the same simulated
#   responses (because all the samples from GA RNAseq are already included, we
#   don't have to simulate all new data sets for CNV and proteomics). Therefore,
#   we can also keep the AES-PCA results for the CNV and proteomics data sets
#   for the MiniMax approach. Question:
#   1. should we use matched-sample techniques AND feature-matched (stacked)
#      pathwayPCA for comparison, or
#   2. should we just use matched-sample techniques for comparison?
#   This could realistically be two different simulation studies, 1) using 222
#   shared samples (which lets us compare to MOGSA, iProFun, etc), and 2) using
#   5 shared samples but only comparing against stacked pathwayPCA. If we keep
#   them completely separate (rather than attempting to "nest" the studies),
#   then we don't have to match the features.
# NOTE: now that I'm thinking about it, the genes between the GA RNAseq and 
#   HiSeq RNAseq aren't the same and they don't share the same intersection.
#   Recall that when we used the HiSeq RNAseq, we had 5252 genes shared across
#   the three platforms. Now, using GA, we only have 1710. This means that if
#   we want to be able to use the same results for CNV and proteomics, we would
#   have to simulate the missing genes for GA. This could be a problem.

# We will have to simulate fresh data. Also, we only need to save the 1710 genes
#   in the intersection of the three platforms (this will save on space).
