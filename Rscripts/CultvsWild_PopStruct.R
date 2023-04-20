### The BEGINNING ~~~~~
##
# ~ Plots --PCA, admixture,  | First written by Nicolas Lou with later modifications by Homère J. Alves Monteiro
rm(list = ls(all = TRUE))
# Loads functions
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions.R")
# Loads required packages ~
pacman::p_load(tidyverse, cowplot, knitr, scales, MetBrewer,RcppCNPy)

# Reorders Population ~
annot <- read.table("~/Desktop/Scripts/EUostrea/01_infofiles/bamlist_EUostrea.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
annot$V2 <- factor(annot$V2, ordered = T,
                   levels = c("MOLU", "ZECE", "CRES","ORIS","CORS", "PONT",  "RIAE", "MORL",
                              "USAM", "TOLL", "COLN", "BARR","TRAL", "CLEW", "RYAN",
                              "GREV", "WADD", "NISS","LOGS","VENO", "HALS", "THIS",
                              "KALV", "HYPP","LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))




# Inv Reg scaffold4 
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/InvReg_EUostrea/7feb23_scaffold4_InvReg_pcangsd.cov")) 
pca_Reg04 <- PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F)
