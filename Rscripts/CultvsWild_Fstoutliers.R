### The BEGINNING ~~~~~
##
# ~ Plots  | By Nicolas R. Lou and George Pacheco modified by Homère J. Alves Monteiro


# Cleans the environment ~ 
rm(list=ls())
library(tidyverse)
library(cowplot)
library(ggrepel)
library(knitr)
library(RcppCNPy)
library(IRanges)
library(GenomicRanges)
# Loads required packages ~
pacman::p_load(tidyverse, extrafont, lemon, data.table, MetBrewer, cowplot, GenomicRanges, IRanges, RcppCNPy, ggrepel, seqinr)

###### Load and format the data ########
setwd("~/Desktop/Scripts/Data/Cult.vs.Wild/")

# Loads datasets ~ mindInd0.1 15kb window 15kb steps Cult vs Wild dataset
NELL.RYAN_CvsW <- read.table("Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/Sept22_Cult.vs.Wild_mindIND0.1_NELL.RYAN_15KB_15KB--Fst.tsv", header = FALSE)
#NELL.SYDK_CvsW <- read.table("Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/Sept22_Cult.vs.Wild_mindIND0.1_NELL.SYDK_15KB_15KB--Fst.tsv", header = FALSE)
SYDK.RAMS_CvsW <- read.table("Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/Sept22_Cult.vs.Wild_mindIND0.1_SYDK.RAMS_15KB_15KB--Fst.tsv", header = FALSE)
#SYDK.KALV_CvsW <- read.table("Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/Sept22_Cult.vs.Wild_mindIND0.1_SYDK.KALV_15KB_15KB--Fst.tsv", header = FALSE)

colnames(NELL.RYAN_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(SYDK.RAMS_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
#colnames(NELL.SYDK_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
#colnames(SYDK.KALV_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")

NELL.RYAN_CvsW <- NELL.RYAN_CvsW[order(as.numeric(substr(NELL.RYAN_CvsW$CHR, 9, nchar(NELL.RYAN_CvsW$CHR)))), ] #super important
NELL.RYAN_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs. Loch Ryan W (SCO)"))
NELL.RYAN_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.RYAN_CvsW))

SYDK.RAMS_CvsW <- SYDK.RAMS_CvsW[order(as.numeric(substr(SYDK.RAMS_CvsW$CHR, 9, nchar(SYDK.RAMS_CvsW$CHR)))), ] #super important
SYDK.RAMS_CvsW$Species <- factor(paste("Sydkoster C  (SWE) vs. Ramsholmen W (SWE)"))
SYDK.RAMS_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.RAMS_CvsW))

# #NELL.SYDK_CvsW <- NELL.SYDK_CvsW[order(as.numeric(substr(NELL.SYDK_CvsW$CHR, 9, nchar(NELL.SYDK_CvsW$CHR)))), ] #super important
# #NELL.SYDK_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs. Sydkoster C (SWE)"))
# NELL.SYDK_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.SYDK_CvsW))
# 
# SYDK.KALV_CvsW <- SYDK.KALV_CvsW[order(as.numeric(substr(SYDK.KALV_CvsW$CHR, 9, nchar(SYDK.KALV_CvsW$CHR)))), ] #super important
# SYDK.KALV_CvsW$Species <- factor(paste(" Sydkoster C (SWE) vs. Kalvö W (SWE) "))
# SYDK.KALV_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.KALV_CvsW))

NELL.RYAN_CvsW_CHR_IDs <- as.data.frame(unique(NELL.RYAN_CvsW$CHR)); colnames(NELL.RYAN_CvsW_CHR_IDs) <- c("CHR")
NELL.RYAN_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.RYAN_CvsW_CHR_IDs))
NELL.RYAN_CvsWUp <- merge(NELL.RYAN_CvsW, NELL.RYAN_CvsW_CHR_IDs, by = "CHR")
NELL.RYAN_CvsWUp <- NELL.RYAN_CvsWUp %>% arrange(CHR_IDs)
NELL.RYAN_CvsWUp$CHR_State <- ifelse(NELL.RYAN_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
SYDK.RAMS_CvsW_CHR_IDs <- as.data.frame(unique(SYDK.RAMS_CvsW$CHR)); colnames(SYDK.RAMS_CvsW_CHR_IDs) <- c("CHR")
SYDK.RAMS_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(SYDK.RAMS_CvsW_CHR_IDs))
SYDK.RAMS_CvsWUp <- merge(SYDK.RAMS_CvsW, SYDK.RAMS_CvsW_CHR_IDs, by = "CHR")
SYDK.RAMS_CvsWUp <- SYDK.RAMS_CvsWUp %>% arrange(CHR_IDs)
SYDK.RAMS_CvsWUp$CHR_State <- ifelse(SYDK.RAMS_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
# NELL.SYDK_CvsW_CHR_IDs <- as.data.frame(unique(NELL.SYDK_CvsW$CHR)); colnames(NELL.SYDK_CvsW_CHR_IDs) <- c("CHR")
# NELL.SYDK_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.SYDK_CvsW_CHR_IDs))
# NELL.SYDK_CvsWUp <- merge(NELL.SYDK_CvsW, NELL.SYDK_CvsW_CHR_IDs, by = "CHR")
# NELL.SYDK_CvsWUp <- NELL.SYDK_CvsWUp %>% arrange(CHR_IDs)
# NELL.SYDK_CvsWUp$CHR_State <- ifelse(NELL.SYDK_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
# SYDK.KALV_CvsW_CHR_IDs <- as.data.frame(unique(SYDK.KALV_CvsW$CHR)); colnames(SYDK.KALV_CvsW_CHR_IDs) <- c("CHR")
# SYDK.KALV_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(SYDK.KALV_CvsW_CHR_IDs))
# SYDK.KALV_CvsWUp <- merge(SYDK.KALV_CvsW, SYDK.KALV_CvsW_CHR_IDs, by = "CHR")
# SYDK.KALV_CvsWUp <- SYDK.KALV_CvsWUp %>% arrange(CHR_IDs)
# SYDK.KALV_CvsWUp$CHR_State <- ifelse(SYDK.KALV_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")

# Gets column names ~
#fulldf <- rbind(NELL.RYANUp,CLEW.NELLUp, CLEW.RYANUp)
#fulldf <- rbind(NELL.RYAN_CvsWUp, NELL.SYDK_CvsWUp, SYDK.RAMS_CvsWUp, SYDK.KALV_CvsWUp)
fulldf <- rbind(NELL.RYAN_CvsWUp, SYDK.RAMS_CvsWUp)


###### Load and format the data 2 ########
setwd("~/Desktop/Scripts/Data/Cult.vs.Wild/")

# Loads datasets ~ SLFst5kb/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23 Cult vs Wild dataset
NELL.RYAN_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_NELL.RYAN_5KB_5KB--Fst.tsv", header = FALSE)
# NELL.TRAL_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_NELL.TRAL_5KB_5KB--Fst.tsv", header = FALSE)
# NELL.SYDK_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_NELL.SYDK_5KB_5KB--Fst.tsv", header = FALSE)
# NELL.RAMS_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_NELL.RAMS_5KB_5KB--Fst.tsv", header = FALSE)
# NELL.KALV_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_NELL.KALV_5KB_5KB--Fst.tsv", header = FALSE)
SYDK.RAMS_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_SYDK.RAMS_5KB_5KB--Fst.tsv", header = FALSE)
# SYDK.KALV_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_SYDK.KALV_5KB_5KB--Fst.tsv", header = FALSE)
# SYDK.RYAN_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_RYAN.SYDK_5KB_5KB--Fst.tsv", header = FALSE)
# SYDK.TRAL_CvsW <- read.table("SLFst5kb_propersetMinDepth100Jan23/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_SYDK.TRAL_5KB_5KB--Fst.tsv", header = FALSE)



colnames(NELL.RYAN_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
# colnames(NELL.TRAL_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
# colnames(NELL.SYDK_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
# colnames(NELL.RAMS_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
# colnames(NELL.KALV_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(SYDK.RAMS_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
# colnames(SYDK.KALV_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
# colnames(SYDK.RYAN_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
# colnames(SYDK.TRAL_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")

NELL.RYAN_CvsW <- NELL.RYAN_CvsW[order(as.numeric(substr(NELL.RYAN_CvsW$CHR, 9, nchar(NELL.RYAN_CvsW$CHR)))), ] #super important
NELL.RYAN_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs. Loch Ryan W (SCO)"))
NELL.RYAN_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.RYAN_CvsW))

# NELL.TRAL_CvsW <- NELL.TRAL_CvsW[order(as.numeric(substr(NELL.TRAL_CvsW$CHR, 9, nchar(NELL.TRAL_CvsW$CHR)))), ] #super important
# NELL.TRAL_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs. Tralee Bay W (IRE)"))
# NELL.TRAL_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.TRAL_CvsW))

# NELL.SYDK_CvsW <- NELL.SYDK_CvsW[order(as.numeric(substr(NELL.SYDK_CvsW$CHR, 9, nchar(NELL.SYDK_CvsW$CHR)))), ] #super important
# NELL.SYDK_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs.  Sydkoster C (SWE)"))
# NELL.SYDK_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.SYDK_CvsW))
# 
# NELL.RAMS_CvsW <- NELL.RAMS_CvsW[order(as.numeric(substr(NELL.RAMS_CvsW$CHR, 9, nchar(NELL.RAMS_CvsW$CHR)))), ] #super important
# NELL.RAMS_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs.  Ramsholmen W (SWE)"))
# NELL.RAMS_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.RAMS_CvsW))
# 
# NELL.KALV_CvsW <- NELL.KALV_CvsW[order(as.numeric(substr(NELL.KALV_CvsW$CHR, 9, nchar(NELL.KALV_CvsW$CHR)))), ] #super important
# NELL.KALV_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs.  Kalvö W (SWE)"))
# NELL.KALV_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.KALV_CvsW))

SYDK.RAMS_CvsW <- SYDK.RAMS_CvsW[order(as.numeric(substr(SYDK.RAMS_CvsW$CHR, 9, nchar(SYDK.RAMS_CvsW$CHR)))), ] #super important
SYDK.RAMS_CvsW$Species <- factor(paste("Sydkoster C  (SWE) vs. Ramsholmen W (SWE)"))
SYDK.RAMS_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.RAMS_CvsW))

# SYDK.KALV_CvsW <- SYDK.KALV_CvsW[order(as.numeric(substr(SYDK.KALV_CvsW$CHR, 9, nchar(SYDK.KALV_CvsW$CHR)))), ] #super important
# SYDK.KALV_CvsW$Species <- factor(paste(" Sydkoster C (SWE) vs. Kalvö W (SWE) "))
# SYDK.KALV_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.KALV_CvsW))
# 
# SYDK.RYAN_CvsW <- SYDK.RYAN_CvsW[order(as.numeric(substr(SYDK.RYAN_CvsW$CHR, 9, nchar(SYDK.RYAN_CvsW$CHR)))), ] #super important
# SYDK.RYAN_CvsW$Species <- factor(paste(" Sydkoster C (SWE) vs. Loch Ryan W (SCO)"))
# SYDK.RYAN_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.RYAN_CvsW))
# 
# SYDK.TRAL_CvsW <- SYDK.TRAL_CvsW[order(as.numeric(substr(SYDK.TRAL_CvsW$CHR, 9, nchar(SYDK.TRAL_CvsW$CHR)))), ] #super important
# SYDK.TRAL_CvsW$Species <- factor(paste(" Sydkoster C (SWE) vs. Tralee Bay W (IRE)"))
# SYDK.TRAL_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.TRAL_CvsW))


NELL.RYAN_CvsW_CHR_IDs <- as.data.frame(unique(NELL.RYAN_CvsW$CHR)); colnames(NELL.RYAN_CvsW_CHR_IDs) <- c("CHR")
NELL.RYAN_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.RYAN_CvsW_CHR_IDs))
NELL.RYAN_CvsWUp <- merge(NELL.RYAN_CvsW, NELL.RYAN_CvsW_CHR_IDs, by = "CHR")
NELL.RYAN_CvsWUp <- NELL.RYAN_CvsWUp %>% arrange(CHR_IDs)
NELL.RYAN_CvsWUp$CHR_State <- ifelse(NELL.RYAN_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")

# NELL.TRAL_CvsW
# NELL.TRAL_CvsW_CHR_IDs <- as.data.frame(unique(NELL.TRAL_CvsW$CHR)); colnames(NELL.TRAL_CvsW_CHR_IDs) <- c("CHR")
# NELL.TRAL_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.TRAL_CvsW_CHR_IDs))
# NELL.TRAL_CvsWUp <- merge(NELL.TRAL_CvsW, NELL.TRAL_CvsW_CHR_IDs, by = "CHR")
# NELL.TRAL_CvsWUp <- NELL.TRAL_CvsWUp %>% arrange(CHR_IDs)
# NELL.TRAL_CvsWUp$CHR_State <- ifelse(NELL.TRAL_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
# 
# NELL.SYDK_CvsW_CHR_IDs <- as.data.frame(unique(NELL.SYDK_CvsW$CHR)); colnames(NELL.SYDK_CvsW_CHR_IDs) <- c("CHR")
# NELL.SYDK_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.SYDK_CvsW_CHR_IDs))
# NELL.SYDK_CvsWUp <- merge(NELL.SYDK_CvsW, NELL.SYDK_CvsW_CHR_IDs, by = "CHR")
# NELL.SYDK_CvsWUp <- NELL.SYDK_CvsWUp %>% arrange(CHR_IDs)
# NELL.SYDK_CvsWUp$CHR_State <- ifelse(NELL.SYDK_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
# 
# NELL.RAMS_CvsW_CHR_IDs <- as.data.frame(unique(NELL.RAMS_CvsW$CHR)); colnames(NELL.RAMS_CvsW_CHR_IDs) <- c("CHR")
# NELL.RAMS_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.RAMS_CvsW_CHR_IDs))
# NELL.RAMS_CvsWUp <- merge(NELL.RAMS_CvsW, NELL.RAMS_CvsW_CHR_IDs, by = "CHR")
# NELL.RAMS_CvsWUp <- NELL.RAMS_CvsWUp %>% arrange(CHR_IDs)
# NELL.RAMS_CvsWUp$CHR_State <- ifelse(NELL.RAMS_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
# 
# NELL.KALV_CvsW
# NELL.KALV_CvsW_CHR_IDs <- as.data.frame(unique(NELL.KALV_CvsW$CHR)); colnames(NELL.KALV_CvsW_CHR_IDs) <- c("CHR")
# NELL.KALV_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.RAMS_CvsW_CHR_IDs))
# NELL.KALV_CvsWUp <- merge(NELL.KALV_CvsW, NELL.KALV_CvsW_CHR_IDs, by = "CHR")
# NELL.KALV_CvsWUp <- NELL.KALV_CvsWUp %>% arrange(CHR_IDs)
# NELL.KALV_CvsWUp$CHR_State <- ifelse(NELL.KALV_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")

SYDK.RAMS_CvsW_CHR_IDs <- as.data.frame(unique(SYDK.RAMS_CvsW$CHR)); colnames(SYDK.RAMS_CvsW_CHR_IDs) <- c("CHR")
SYDK.RAMS_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(SYDK.RAMS_CvsW_CHR_IDs))
SYDK.RAMS_CvsWUp <- merge(SYDK.RAMS_CvsW, SYDK.RAMS_CvsW_CHR_IDs, by = "CHR")
SYDK.RAMS_CvsWUp <- SYDK.RAMS_CvsWUp %>% arrange(CHR_IDs)
SYDK.RAMS_CvsWUp$CHR_State <- ifelse(SYDK.RAMS_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")

# SYDK.KALV_CvsW_CHR_IDs <- as.data.frame(unique(SYDK.KALV_CvsW$CHR)); colnames(SYDK.KALV_CvsW_CHR_IDs) <- c("CHR")
# SYDK.KALV_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(SYDK.KALV_CvsW_CHR_IDs))
# SYDK.KALV_CvsWUp <- merge(SYDK.KALV_CvsW, SYDK.KALV_CvsW_CHR_IDs, by = "CHR")
# SYDK.KALV_CvsWUp <- SYDK.KALV_CvsWUp %>% arrange(CHR_IDs)
# SYDK.KALV_CvsWUp$CHR_State <- ifelse(SYDK.KALV_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
# 
# SYDK.RYAN_CvsW_CHR_IDs
# SYDK.RYAN_CvsW_CHR_IDs <- as.data.frame(unique(SYDK.RYAN_CvsW$CHR)); colnames(SYDK.RYAN_CvsW_CHR_IDs) <- c("CHR")
# SYDK.RYAN_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(SYDK.RYAN_CvsW_CHR_IDs))
# SYDK.RYAN_CvsWUp <- merge(SYDK.RYAN_CvsW, SYDK.RYAN_CvsW_CHR_IDs, by = "CHR")
# SYDK.RYAN_CvsWUp <- SYDK.RYAN_CvsWUp %>% arrange(CHR_IDs)
# SYDK.RYAN_CvsWUp$CHR_State <- ifelse(SYDK.RYAN_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
# 
# SYDK.TRAL_CvsW_CHR_IDs <- as.data.frame(unique(SYDK.TRAL_CvsW$CHR)); colnames(SYDK.TRAL_CvsW_CHR_IDs) <- c("CHR")
# SYDK.TRAL_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(SYDK.TRAL_CvsW_CHR_IDs))
# SYDK.TRAL_CvsWUp <- merge(SYDK.TRAL_CvsW, SYDK.TRAL_CvsW_CHR_IDs, by = "CHR")
# SYDK.TRAL_CvsWUp <- SYDK.TRAL_CvsWUp %>% arrange(CHR_IDs)
# SYDK.TRAL_CvsWUp$CHR_State <- ifelse(SYDK.TRAL_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")

# Gets column names ~
#fulldf <- rbind(NELL.RYANUp,CLEW.NELLUp, CLEW.RYANUp)
# fulldf_sco <- rbind(NELL.RYAN_CvsWUp, NELL.TRAL_CvsWUp, NELL.RAMS_CvsWUp, NELL.KALV_CvsWUp, NELL.SYDK_CvsWUp)
# fulldf_swe <- rbind(NELL.SYDK_CvsWUp, SYDK.RAMS_CvsWUp, SYDK.KALV_CvsWUp, SYDK.RYAN_CvsWUp, SYDK.TRAL_CvsWUp)
fulldf_5kb <- rbind(NELL.RYAN_CvsWUp, SYDK.RAMS_CvsWUp)



###### Create plots ########
Fst_Plot_CultvsWild_20oct22 <-
  ggplot() +
  geom_point(data = fulldf, aes(x = gPoint_c, y = Fst, fill= CHR_State, colour = CHR_State), shape = 21, size = 1, alpha = 0.4) +
  facet_rep_grid(Species~. , scales = "free_x") +
  scale_x_continuous("Chromosomes",
                     expand = c(.005, .005)) +
  scale_y_continuous("Fst (5Kb Sliding Windows)",
                     breaks = c(.30, .60, .90), 
                     labels = c(".30", ".60", ".90"),
                     limits = c(0, .99),
                     expand = c(0.01, 0.01)) +
  scale_fill_manual(values = MetBrewer::met.brewer("OKeeffe1", n = 2, type = "discrete")) +
  scale_colour_manual(values = MetBrewer::met.brewer("OKeeffe1", n = 2, type = "discrete")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text.x = element_text(colour = "#000000", size = 6),
        axis.text.y = element_text(colour = "#000000", size = 12, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", size = .3),
        axis.ticks.y = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .5),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = "none", fill = "none")
last_plot()

ggsave(Fst_Plot_CultvsWild_20oct22, file = "~/Desktop/Scripts/Cult.vs.Wild/Figures/Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/RYANvsNELL_SYDKvsRAMS_15kb15kb_mindind0.1.pdf",
       device = cairo_pdf, scale = 1, width = 30, height = 15, dpi = 300)
dev.off()


# # Calculate the top 5% of points based on their Fst values
# top1pct <- fulldf %>% 
#   arrange(desc(Fst)) %>% 
#   slice(1:round(nrow(.) * 0.01))
# # calculate the top 0.5% of points based on their Fst values by species separately
# top0.1pct_species <- fulldf %>% 
#   group_by(Species) %>% 
#   arrange(desc(Fst)) %>% 
#   slice(1:round(nrow(.) * 0.001))
# 
# # Create the plot
# top0.1pct_speciesplot <- ggplot() +
#   geom_point(data = fulldf, aes(x = gPoint_c, y = Fst, fill= CHR_State, colour = CHR_State), shape = 21, size = .1, alpha = 0.4) +
#   geom_point(data = top0.1pct_species, aes(x = gPoint_c, y = Fst), shape = 21, size = 2, alpha = 1) +
#   facet_rep_grid(Species~. , scales = "free_x") +
#   scale_x_continuous("Chromosomes",
#                      expand = c(.005, .005)) +
#   scale_y_continuous("Fst (15Kb Sliding Windows)",
#                      breaks = c(.30, .60, .90), 
#                      labels = c(".30", ".60", ".90"),
#                      limits = c(0, .99),
#                      expand = c(0.01, 0.01)) +
#   scale_colour_manual(values = MetBrewer::met.brewer("OKeeffe1", n = 2, type = "discrete"))+
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "#000000", size = .3),
#         axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
#         axis.text.x = element_text(colour = "#000000", size = 6),
#         axis.text.y = element_text(colour = "#000000", size = 12, face = "bold"),
#         axis.ticks.x = element_line(color = "#000000", size = .3),
#         axis.ticks.y = element_line(color = "#000000", size = .3),
#         strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .5),
#         strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
#         legend.position = "top",
#         legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
#         legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
#         legend.key = element_rect(fill = NA),
#         legend.background = element_blank()) +
#   guides(colour = "none", fill = "none")
# last_plot()
# ggsave(top0.1pct_speciesplot, file = "~/Desktop/Scripts/Cult.vs.Wild/Figures/Top_1percent_Fst_perComparison_15kbwin_15kbsteps_CultvsWildpops_minind0.1RYANvsNELL_SYDKvsRAMS_15kb15kb_mindind0.1.pdf",
#        device = cairo_pdf, scale = 1, width = 30, height = 15, dpi = 300)
# dev.off()


# Create the ggplot object and add the data
 
top01pctplot <- ggplot() +
  geom_point(data = fulldf, aes(x = gPoint_c, y = Fst, fill = CHR_State, colour = CHR_State), shape = 21, size = .1, alpha = 0.4) +
  
  # Calculate the top 1% of Fst for each species and add the points to the plot
  geom_point(data = fulldf %>% group_by(Species) %>% top_n(n = floor(n() * 0.001), Fst),
             aes(x = gPoint_c, y = Fst), shape = 21, size = 2, alpha = 1) +
  # Create a grid of facets for each species
  facet_rep_grid(Species~., scales = "free_x") +
  # Add labels and customize the scales
  scale_x_continuous("Chromosomes", expand = c(.005, .005)) +
  scale_y_continuous("Fst (15Kb Sliding Windows)",
                     breaks = c(.30, .60, .90),
                     labels = c(".30", ".60", ".90"),
                     limits = c(0, .99),
                     expand = c(0.01, 0.01)) +
  scale_colour_manual(values = MetBrewer::met.brewer("OKeeffe1", n = 2, type = "discrete"))+
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text.x = element_text(colour = "#000000", size = 6),
        axis.text.y = element_text(colour = "#000000", size = 12, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", size = .3),
        axis.ticks.y = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .4),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = "none", fill = "none")
last_plot()
ggsave(top01pctplot, file = "~/Desktop/Scripts/Cult.vs.Wild/Figures/Top_0.1percent_Fst_15kbwin15kbsteps_CultvsWildmindind0.1.pdf",
       device = cairo_pdf, scale = 1, width = 30, height = 15, dpi = 300)
dev.off()


# Find and output some of the most obvious Fst outliers
fulldf_outliers <- fulldf %>%
  filter(CHR_IDs %in% seq(1:10)) %>%
  group_by(CHR_IDs) %>%
  mutate(outlier_status = ifelse(Fst > 0.3, T, F))
fulldf_outliers_ID <- fulldf_outliers %>%
  filter(outlier_status == T) %>%
  summarize(min_pos = min(SNP), max_pos = max(SNP))
fulldf_outliers_ID
# Find and output some of the most obvious Fst outliers with Fst > 0.6 for LOch nell comparisons and Fst > 0.3 for Sydkoster comparisons.
fulldf_outliers <- fulldf %>%
  filter(CHR_IDs %in% seq(1:10)) %>%
  group_by(CHR_IDs) %>%
  mutate(outlier_status = ifelse(Fst > 0.3 & !(Species %in% c("Loch Nell C (SCO) vs. Loch Ryan W (SCO)", "Loch Nell C (SCO) vs. Sydkoster C (SWE)")), T,
                                 ifelse(Fst > 0.6 & (Species %in% c("Loch Nell C (SCO) vs. Loch Ryan W (SCO)", "Loch Nell C (SCO) vs. Sydkoster C (SWE)")), T,
                                        ifelse(Fst > 0.3 & (Species %in% c("Sydkoster C (SWE) vs. Ramsholmen W (SWE)", "Sydkoster C (SWE) vs. Kalvö W (SWE)")), T, F))))
  
# List of the outliers
fulldf_outliers_ID <- fulldf_outliers %>%
  filter(outlier_status == T)
fulldf_outliers_ID%>% kable()
write_tsv(fulldf_outliers_ID, "~/Desktop/Scripts/Cult.vs.Wild/Figures/17Apr23--LIST_obviousFstoutliers.tsv")


#Create a new variable that is a scaled version of gPoint_c
fulldf_outliers$scaled_gPoint_c <- as.numeric(as.factor(fulldf_outliers$gPoint_c))
#Plot outliers on chr 1,2 and 10 using the scaled variabl
Outlier_plot2 <- ggplot(fulldf_outliers, aes(x=scaled_gPoint_c, y=Fst, color=outlier_status, shape=CHR_IDs)) +
  geom_point(size=ifelse(fulldf_outliers$outlier_status == "FALSE", 1, 4),
             color=ifelse(fulldf_outliers$outlier_status == "FALSE", "black", "orange"),
             aes(shape = CHR)) +
  scale_x_discrete("Pos",
                   expand = c(.005, .005)) +
  facet_wrap(~Species, ncol = 1) +  
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))+
  theme_cowplot()+
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text.x = element_text(colour = "#000000", size = 6),
        axis.text.y = element_text(colour = "#000000", size = 12, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", size = .3),
        axis.ticks.y = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .4),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = "none", fill = "none")

Outlier_plot3 <-ggplot(fulldf_outliers, aes(x=scaled_gPoint_c, y=Fst, color=outlier_status, shape=CHR_IDs)) +
  geom_point(size=ifelse(fulldf_outliers$outlier_status == "FALSE", 1, 4),
             color=ifelse(fulldf_outliers$outlier_status == "FALSE", "black", "orange"),
             aes(shape = CHR),
             position = position_dodge(width = 0.3)) +
  scale_x_discrete("Pos",
                   expand = c(.005, .005)) +
  facet_wrap(~Species, ncol = 1) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  theme_cowplot() +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text.x = element_text(colour = "#000000", size = 6),
        axis.text.y = element_text(colour = "#000000", size = 12, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", size = .3),
        axis.ticks.y = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .4),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = "none", fill = "none")



ggsave(Outlier_plot2, file = "~/Desktop/Scripts/Cult.vs.Wild/Figures/Apr23--obviousFstoutliers_FstMan.pdf",
       device = cairo_pdf, scale = 1, width = 30, height = 15, dpi = 300)
ggsave(Outlier_plot2, file = "~/Desktop/Scripts/Cult.vs.Wild/Figures/Apr23--obviousFstoutliers_FstMan.png"
       , scale = 1, width = 30, height = 15, dpi = 300)
ggsave(Outlier_plot3, file = "~/Desktop/Scripts/Cult.vs.Wild/Figures/Apr23_alt--obviousFstoutliers_FstMan.png"
       , scale = 1, width = 30, height = 15, dpi = 300)


dev.off()


#ggplot(fulldf_outliers, aes(Fst)) + geom_histogram(bins=100)

fulldf_outliers_ID <- fulldf_outliers %>%
  filter(outlier_status == T) %>%
  dplyr::summarize(min_pos = min(SNP), max_pos = max(SNP))
## note that more outliers can be added to this list later
fulldf_outliers_ID %>% kable()
fulldf_outliers_ID
#Plot the distribution of Fst values
ggplot(fulldf, aes(Fst)) + geom_histogram(bins=100)
# Print the dataframe after mutating outlier_status
fulldf %>% filter(CHR_IDs %in% seq(1:10)) %>% 
    group_by(CHR_IDs) %>%
    mutate(outlier_status = ifelse(Fst > 0.3, T, F))


#write the outliers
## Extract these outliers and output them as fasta files
Oedu <- read_lines("../../../IGV/fileOegenome10scaffoldC3G.fasta")
lg01_sequence <- Oedu[(grep("scaffold1", Oedu)-1):grep("scaffold2", Oedu)] %>% str_c(collapse = "")
lg01_outlier <- str_sub(lg01_sequence, 89242500, 106387500)

lg02_sequence <- Oedu[(grep("scaffold2", Oedu)-1):grep("scaffold3", Oedu)] %>% str_c(collapse = "")
lg02_outlier <- str_sub(lg02_sequence, 24442500, 24452500)

lg10_sequence <-  str_extract(Oedu, "(?<=scaffold10)[^\n]+") %>% str_trim()
lg10_sub_sequence <- sub(".*scaffold10(.*?)\n.*", "\\1", lg10_sequence)
lg10_sub_sequence <- lg10_sub_sequence[!is.na(lg10_sub_sequence)]


# read the fasta file
Oedu2 <- read.fasta("../../../IGV/fileOegenome10scaffoldC3G.fasta")
# extract the scaffold10 sequence
scaffold10_sequence <- Oedu[grep("scaffold10", names(Oedu))]
lg10_outlier <- str_sub(lg10_sequence, 11572500, 28162500)


outlier_fasta_chr2 <- c(
  ">LG02_outlier",
  lg02_outlier)

write_lines(outlier_fasta_chr2, "../Cult.vs.Wild/chr2_NELLvsRYANSYDK_fst_outlier.fasta")

###### Functional enrichment analysis ###### 
####Functions####
#Define a function to find genes close to outlier SNPs
find_nearby_genes <- function(outliers, annotation, distance){
  annotation_range <- with(annotation, GRanges(lg, IRanges(pos_start, pos_end), strand))
  for (i in seq_len(dim(outliers)[1])){
    outlier <- outliers[i,]
    outlier_lg_pos_bin <- outlier$lg_pos_bin
    outlier_range <- with(outlier, GRanges(lg, IRanges(min_pos - distance, max_pos + distance)))
    overlapping_genes_index <- countOverlaps(annotation_range, outlier_range)>0
    overlapping_genes <- annotation[overlapping_genes_index, ] %>%
      mutate(lg_pos_bin = outlier_lg_pos_bin)
    if (i==1){
      overlapping_genes_final <- overlapping_genes
    } else {
      overlapping_genes_final <- bind_rows(overlapping_genes_final, overlapping_genes)
    }
  }
  return(overlapping_genes_final)
}
#Define a function to plot outlier position and nearby genes
plot_nearby_genes <- function(outlier_sites, nearby_genes){
  ggplot(nearby_genes, aes(ymin=pos_start, ymax=pos_end, y=(pos_start+pos_end)/2, x=0)) +
    geom_linerange(size=1, position=position_jitter(1, 0, 1)) +
    geom_label(aes(label=id, x=-0.2), size=5, position=position_jitter(1, 0, 1)) +
    geom_linerange(data=outlier_sites, aes(ymin=min_pos, ymax=max_pos, x=1, y=(min_pos+max_pos)/2), color="red", size=5) +
    facet_wrap(~lg_pos_bin, ncol = 1,scales="free") +
    coord_flip() +
    theme_cowplot()
}

plot_nearby_genes <- function(outlier_sites, nearby_genes){
  lg_labeller <- function(lg_pos_bin){
    return(str_extract(lg_pos_bin, "[^(]+"))
  }
  mutate(nearby_genes, abbr=ifelse(nchar(abbr>10), str_extract(abbr, "[^ ]+"), abbr), abbr=ifelse(abbr=="unknown", NA, abbr)) %>%
    ggplot() +
    geom_segment(aes(x=pos_start, xend=pos_end, y=strand, yend=strand), size=1, arrow = arrow(length = unit(0.2,"cm"))) +
    geom_segment(aes(x=pos_end, xend=pos_start, y=strand, yend=strand), size=1, arrow = arrow(length = unit(0.2,"cm"))) +
    geom_text_repel(aes(label=id, x=(pos_start+pos_end)/2, y=strand), size=4, position = position_nudge(y=0.5)) +
    geom_segment(data=outlier_sites, aes(x=min_pos, xend=max_pos, y="outlier", yend="outlier"), color="red", size=2, arrow = arrow(length = unit(0.2,"cm"))) +
    geom_segment(data=outlier_sites, aes(x=max_pos, xend=min_pos, y="outlier", yend="outlier"), color="red", size=2, arrow = arrow(length = unit(0.2,"cm"))) +
    facet_wrap(~lg_pos_bin, ncol = 1, scales="free_x", strip.position="right", labeller = labeller(lg_pos_bin = lg_labeller)) +
    xlab("position") +
    theme_cowplot()+
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "#000000", size = .3),
          axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
          axis.text.x = element_text(colour = "#000000", size = 6),
          axis.text.y = element_text(colour = "#000000", size = 12, face = "bold"),
          axis.ticks.x = element_line(color = "#000000", size = .3),
          axis.ticks.y = element_line(color = "#000000", size = .3),
          strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .4),
          strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
          legend.position = "top",
          legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
          legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
          legend.key = element_rect(fill = NA),
          legend.background = element_blank()) +
    guides(colour = "none", fill = "none")
  
    
}



## Only run this once
annotation_gene_MOR <- read_tsv("/Users/homere/Desktop/IGV/gff_oedulis.tsv", skip = 1, col_names = F) %>%
  filter(X3=="gene")

#annotation_gene_ROS <- read_tsv("/Users/homere/Desktop/IGV/GCF_023158985.2_OEROSLIN.tsv.gff", skip = 1, col_names = F) %>%
#  filter(X3=="gene")


write_tsv(annotation_gene_MOR, "/Users/homere/Desktop/IGV/Oegenome10chromassembly_MOR_gene_only.tsv")
# write_tsv(annotation_gene_ROS, "/Users/homere/Desktop/IGV/Oegenome10chromassembly_ROS_gene_only.tsv")

annotation_gene_cleaned_MOR <-  read_tsv("/Users/homere/Desktop/IGV/Oegenome10chromassembly_MOR_gene_only.tsv", col_names = T) %>%
  separate(X9, c("id", "name", "alias", "note", "dbxref", "ontology_term", "junk"), sep = ";") %>%
  select(lg=X1, pos_start=X4, pos_end=X5, strand=X7, note=note, ontology_term=ontology_term, id) %>%
  mutate(note=str_remove(note, "Note=")) %>%
  mutate(ontology_term=str_remove(ontology_term, "Ontology_term=")) %>%
  mutate(abbr=str_remove(note, "Similar to ")) %>%
  mutate(abbr=str_remove(abbr, ":.*$")) %>%
  mutate(abbr=ifelse(abbr=="Protein of unknown function", "unknown", abbr))
write_tsv(annotation_gene_cleaned_MOR, "/Users/homere/Desktop/IGV/MOR_gene_only_cleaned.tsv")

# #annotation_gene_cleaned_ROS <-  read_tsv("/Users/homere/Desktop/IGV/Oegenome10chromassembly_ROS_gene_only.tsv", col_names = T) %>%
#   separate(X9, c("id", "name", "alias", "note", "dbxref", "ontology_term", "junk"), sep = ";") %>%
#   select(lg=X1, pos_start=X4, pos_end=X5, strand=X7, note=note, ontology_term=ontology_term) %>%
#   mutate(note=str_remove(note, "Note=")) %>%
#   mutate(ontology_term=str_remove(ontology_term, "Ontology_term=")) %>%
#   mutate(abbr=str_remove(note, "Similar to ")) %>%
#   mutate(abbr=str_remove(abbr, ":.*$")) %>%
#   mutate(abbr=ifelse(abbr=="Protein of unknown function", "unknown", abbr))
# write_tsv(annotation_gene_cleaned_ROS, "/Users/homere/Desktop/IGV/ROS_gene_only_cleaned.tsv")
# unique(annotation_gene_cleaned_ROS$lg)

## Read in the cleaned annotation file
annotation_gene_cleaned_MOR_lg <- read_tsv("/Users/homere/Desktop/IGV/MOR_gene_only_cleaned.tsv") %>%
  filter(str_detect(lg, "scaffold"))
# annotation_gene_cleaned_ROS_lg <- read_tsv("/Users/homere/Desktop/IGV/ROS_gene_only_cleaned.tsv") %>%
#   filter(str_detect(lg, "scaffold"))

#Extract outlier positions
## Read in outliers from Fst outliers manually
# Create a dataframe with the desired values
genome_selection_sites_outliers <- data.frame(lg = c("scaffold1", "scaffold2", "scaffold5",),
                                              pos_bin = c("(10e+06]", "(10e+06]", "(10e+06]"),
                                              min_pos = c(72547500, 25837500, 40312500),
                                              max_pos = c(72592500, 25822500, 40297500),
                                              length = c(15000, 15000, 15000))

#Add a column with a string combining the chromosome and bin names
genome_selection_sites_outliers$lg_pos_bin <- paste0(genome_selection_sites_outliers$lg, genome_selection_sites_outliers$pos_bin)
set.seed(1)


selection_sites_outliers_nearby_genes_200kb <- find_nearby_genes(genome_selection_sites_outliers, annotation_gene_cleaned_MOR_lg, 200000)
set.seed(1)
plot_nearby_genes(genome_selection_sites_outliers, selection_sites_outliers_nearby_genes_200kb)

#
##
### The END ~~~~~

