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
NELL.SYDK_CvsW <- read.table("Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/Sept22_Cult.vs.Wild_mindIND0.1_NELL.SYDK_15KB_15KB--Fst.tsv", header = FALSE)
SYDK.RAMS_CvsW <- read.table("Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/Sept22_Cult.vs.Wild_mindIND0.1_SYDK.RAMS_15KB_15KB--Fst.tsv", header = FALSE)
SYDK.KALV_CvsW <- read.table("Fst_15kbwin_15kbsteps_CultvsWildpops_minind0.1/Sept22_Cult.vs.Wild_mindIND0.1_SYDK.KALV_15KB_15KB--Fst.tsv", header = FALSE)

colnames(NELL.RYAN_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(SYDK.RAMS_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(NELL.SYDK_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(SYDK.KALV_CvsW) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")

NELL.RYAN_CvsW <- NELL.RYAN_CvsW[order(as.numeric(substr(NELL.RYAN_CvsW$CHR, 9, nchar(NELL.RYAN_CvsW$CHR)))), ] #super important
NELL.RYAN_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs. Loch Ryan W (SCO)"))
NELL.RYAN_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.RYAN_CvsW))

SYDK.RAMS_CvsW <- SYDK.RAMS_CvsW[order(as.numeric(substr(SYDK.RAMS_CvsW$CHR, 9, nchar(SYDK.RAMS_CvsW$CHR)))), ] #super important
SYDK.RAMS_CvsW$Species <- factor(paste("Sydkoster C  (SWE) vs. Ramsholmen W (SWE)"))
SYDK.RAMS_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.RAMS_CvsW))

NELL.SYDK_CvsW <- NELL.SYDK_CvsW[order(as.numeric(substr(NELL.SYDK_CvsW$CHR, 9, nchar(NELL.SYDK_CvsW$CHR)))), ] #super important
NELL.SYDK_CvsW$Species <- factor(paste("Loch Nell C (SCO) vs. Sydkoster C (SWE)"))
NELL.SYDK_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(NELL.SYDK_CvsW))

SYDK.KALV_CvsW <- SYDK.KALV_CvsW[order(as.numeric(substr(SYDK.KALV_CvsW$CHR, 9, nchar(SYDK.KALV_CvsW$CHR)))), ] #super important
SYDK.KALV_CvsW$Species <- factor(paste(" Sydkoster C (SWE) vs. Kalvö W (SWE) "))
SYDK.KALV_CvsW$gPoint_c <- seq(15000, by = 15000, length.out = nrow(SYDK.KALV_CvsW))

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
NELL.SYDK_CvsW_CHR_IDs <- as.data.frame(unique(NELL.SYDK_CvsW$CHR)); colnames(NELL.SYDK_CvsW_CHR_IDs) <- c("CHR")
NELL.SYDK_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(NELL.SYDK_CvsW_CHR_IDs))
NELL.SYDK_CvsWUp <- merge(NELL.SYDK_CvsW, NELL.SYDK_CvsW_CHR_IDs, by = "CHR")
NELL.SYDK_CvsWUp <- NELL.SYDK_CvsWUp %>% arrange(CHR_IDs)
NELL.SYDK_CvsWUp$CHR_State <- ifelse(NELL.SYDK_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")
SYDK.KALV_CvsW_CHR_IDs <- as.data.frame(unique(SYDK.KALV_CvsW$CHR)); colnames(SYDK.KALV_CvsW_CHR_IDs) <- c("CHR")
SYDK.KALV_CvsW_CHR_IDs$CHR_IDs <- seq.int(nrow(SYDK.KALV_CvsW_CHR_IDs))
SYDK.KALV_CvsWUp <- merge(SYDK.KALV_CvsW, SYDK.KALV_CvsW_CHR_IDs, by = "CHR")
SYDK.KALV_CvsWUp <- SYDK.KALV_CvsWUp %>% arrange(CHR_IDs)
SYDK.KALV_CvsWUp$CHR_State <- ifelse(SYDK.KALV_CvsWUp$CHR_IDs %% 2 == 0, "Even", "Odd")

# Gets column names ~
#fulldf <- rbind(NELL.RYANUp,CLEW.NELLUp, CLEW.RYANUp)
fulldf <- rbind(NELL.RYAN_CvsWUp, NELL.SYDK_CvsWUp, SYDK.RAMS_CvsWUp, SYDK.KALV_CvsWUp)



###### Create plots ########
Fst_Plot_CultvsWild_20oct22 <-
  ggplot() +
  geom_point(data = fulldf, aes(x = gPoint_c, y = Fst, fill= CHR_State, colour = CHR_State), shape = 21, size = .1, alpha = 0.4) +
  facet_rep_grid(Species~. , scales = "free_x") +
  scale_x_continuous("Chromosomes",
                     expand = c(.005, .005)) +
  scale_y_continuous("Fst (15Kb Sliding Windows)",
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
  filter(CHR_IDs %in% c("1", "2","10")) %>%
  group_by(CHR_IDs) %>%
  mutate(outlier_status = ifelse(Fst > 0.3, T, F))
fulldf_outliers_ID <- fulldf_outliers %>%
  filter(outlier_status == T) %>%
  summarize(min_pos = min(SNP), max_pos = max(SNP))
fulldf_outliers_ID
# Find and output some of the most obvious Fst outliers with Fst > 0.6 for LOch nell comparisons and Fst > 0.3 for Sydkoster comparisons.
fulldf_outliers <- fulldf %>%
  filter(CHR_IDs %in% c("1", "2","10")) %>%
  group_by(CHR_IDs) %>%
  mutate(outlier_status = ifelse(Fst > 0.3 & !(Species %in% c("Loch Nell C (SCO) vs. Loch Ryan W (SCO)", "Loch Nell C (SCO) vs. Sydkoster C (SWE)")), T,
                                 ifelse(Fst > 0.6 & (Species %in% c("Loch Nell C (SCO) vs. Loch Ryan W (SCO)", "Loch Nell C (SCO) vs. Sydkoster C (SWE)")), T,
                                        ifelse(Fst > 0.3 & (Species %in% c("Sydkoster C (SWE) vs. Ramsholmen W (SWE)", "Sydkoster C (SWE) vs. Kalvö W (SWE)")), T, F))))
  
# List of the outliers
fulldf_outliers_ID <- fulldf_outliers %>%
  filter(outlier_status == T)
fulldf_outliers_ID%>% kable()
write_tsv(fulldf_outliers_ID, "~/Desktop/Scripts/Cult.vs.Wild/Figures/Jan23--LIST_obviousFstoutliers_chr1_2_10_clean.tsv")


#Create a new variable that is a scaled version of gPoint_c
fulldf_outliers$scaled_gPoint_c <- as.numeric(as.factor(fulldf_outliers$gPoint_c))
#Plot outliers on chr 1,2 and 10 using the scaled variabl
Outlier_plot2 <- ggplot(fulldf_outliers, aes(x=scaled_gPoint_c, y=Fst, color=outlier_status, shape=CHR_IDs)) +
  geom_point(size=ifelse(fulldf_outliers$outlier_status == "FALSE", 0.2, 4),aes(shape = CHR)) +
  scale_x_discrete("Pos",
                   expand = c(.005, .005)) +
  facet_wrap(~Species, ncol = 1) +  
  scale_colour_manual(values = MetBrewer::met.brewer("Pillement", n = 2, type = "discrete")) +
  theme_cowplot()+
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

ggsave(Outlier_plot2, file = "~/Desktop/Scripts/Cult.vs.Wild/Figures/Jan23--obviousFstoutliers_chr1_2_10_clean.pdf",
       device = cairo_pdf, scale = 1, width = 30, height = 15, dpi = 300)
dev.off()


#ggplot(fulldf_outliers, aes(Fst)) + geom_histogram(bins=100)

fulldf_outliers_ID <- fulldf_outliers %>%
  filter(outlier_status == T) %>%
  summarize(min_pos = min(SNP), max_pos = max(SNP))
## note that more outliers can be added to this list later
fulldf_outliers_ID %>% kable()
fulldf_outliers_ID
#Plot the distribution of Fst values
ggplot(fulldf, aes(Fst)) + geom_histogram(bins=100)
# Print the dataframe after mutating outlier_status
fulldf %>% filter(CHR_IDs %in% c("1", "2","10")) %>% group_by(CHR_IDs) %>% mutate(outlier_status = ifelse(Fst > 0.3, T, F))


#write the outliers
## Extract these outliers and output them as fasta files
Oedu <- read_lines("../../../IGV/fileOegenome10scaffoldC3G.fasta")
lg01_sequence <- Oedu[(grep("scaffold1", Oedu)-1):grep("scaffold2", Oedu)] %>% str_c(collapse = "")
lg01_outlier <- str_sub(lg01_sequence, 89242500, 106387500)

lg02_sequence <- Oedu[(grep("scaffold2", Oedu)-1):grep("scaffold3", Oedu)] %>% str_c(collapse = "")
lg02_outlier <- str_sub(lg21_sequence, 90967500, 104572500)

# read the fasta file
Oedu2 <- read.fasta("../../../IGV/fileOegenome10scaffoldC3G.fasta")
# extract the scaffold10 sequence
scaffold10_sequence <- Oedu[grep("scaffold10", names(Oedu))]
lg10_outlier <- str_sub(lg10_sequence, 11572500, 28162500)

outlier_fasta <- c(
  ">LG01_outlier",
  lg01_outlier,
  ">LG02_outlier",
  lg02_outlier,
  ">LG10_outlier",
  lg10_outlier
)
write_lines(outlier_fasta, "../angsd/popminind5/bam_list_realigned_mincov_filtered_mindp60_maxdp337_minind33_minq20_popminind5_fst_outlier.fasta")
