#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/Cult.vs.Wild
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N pcangsdApr23
#PBS -e pcangsdApr23.err
#PBS -o pcangsdApr23.out
#PBS -l nodes=1:ppn=30:thinnode
#PBS -l walltime=100:00:00
#PBS -l mem=70gb
#PBS -m n
#PBS -r n

#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20210722
module load angsd/0.940

#pcangsd
module load tools computerome_utils/2.0
module load pcangsd/20220330 

BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
THREADS=10
SNP_LIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'
MINDP=73 # Minimum depth filter
MAXDP=221 # Maximum depth filter
MINQ=20 # Minimum quality filter
MINMAF=0.01 # Minimum minor allele frequency filter
MINMAPQ=20 # Minimum mapping quality (alignment score) filter, default value is 20
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput


## Re-run angsd LD Pruned SNPs list minweight0.2 Global dataset 
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221 \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-setMinDepth $MINDP -setMaxDepth $MAXDP \
-minQ $MINQ -minMapQ $MINMAPQ \
-SNP_pval 1e-6 -minMaf $MINMAF \
-P $THREADS \
$EXTRA_ARG -rmTriallelic 0.05 -trim 0 -baq 1 \
-sites $SNP_LIST \
-rf $LG_LIST

wait

## Create a SNP list to use in downstream analyses
gunzip -c $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.mafs.gz | cut -f 1,2,3,4 | tail -n +2 > $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt

angsd sites index $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt

## Also make it in regions format for downstream analyses
cut -f 1,2 $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt | sed 's/\t/:/g' > $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_regions.txt

## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
cut -f1 $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt | sort | uniq > $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.chrs

wait

#pcangsd for LD pruned SNPS
pcangsd -b $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.beagle.gz \
--selection \
--selection_e 30 \
--pcadapt \
--sites_save \
--snp_weights \
--pcadapt \
--maf_save \
--threads 30 \
--admix \
-e 30 \
-o $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30

#pcangsd for Global VC SNPs
pcangsd -b $OUTPUTFOLDER/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.beagle.gz \
--selection \
--threads 30 \
--sites_save \
--minMaf 0.01 \
-e 30 \
-o $OUTPUTFOLDER/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30


#on c2
R
rm(list=ls(all.names = TRUE))
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)
basedir="/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput"  
genome_selection <- npyLoad(paste0(basedir, "/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30.selection.npy"))  
genome_selection_df <- as.data.frame(genome_selection)
write_tsv(genome_selection_df, paste0(basedir, "/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30.selection.tsv"), col_names = F)


genome_selection <- npyLoad(paste0(basedir, "/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30.selection.npy"))  
genome_selection_sites <- read_tsv(paste0(basedir,"/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt"), 
                                   col_names = c("lg", "pos", "major", "minor")) %>%
  bind_cols(read_table(paste0(basedir,"/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30.sites"), col_names = "keep")) %>% 
  filter(keep==1)%>% 
  dplyr::mutate(PC1=genome_selection[,1],
         PC2=genome_selection[,2],
         PC3=genome_selection[,3],
         PC4=genome_selection[,4],
         PC5=genome_selection[,5],
         PC6=genome_selection[,6],
         PC7=genome_selection[,7],
         PC8=genome_selection[,8],
         PC9=genome_selection[,9],
         PC10=genome_selection[,10]) %>%
  pivot_longer(names_to = "pc", values_to = "chi_squared", cols = 6:15) %>%
  dplyr::select(lg, pos, pc, chi_squared) %>%
  dplyr::mutate(chi_squared = ifelse(pc == "PC1", chi_squared$X1, 
                              ifelse(pc == "PC2", chi_squared$X2,
                                     ifelse(pc == "PC3", chi_squared$X3, 
                                            ifelse(pc == "PC4", chi_squared$X4, 
                                                   ifelse(pc == "PC5", chi_squared$X5,
                                                          ifelse(pc == "PC6", chi_squared$X6,
                                                                 ifelse(pc == "PC7", chi_squared$X7,
                                                                        ifelse(pc == "PC8", chi_squared$X8, 
                                                                               ifelse(pc == "PC9", chi_squared$X9,
                                                                                      ifelse(pc == "PC10", chi_squared$X10,NA))))))))))) %>%
  filter(!is.na(chi_squared)) %>%
  dplyr::mutate(neg_log_p_value=-log(1-pchisq(chi_squared, df=1))) 



#reorder chr
genome_selection_sites$lg <- factor(genome_selection_sites$lg , ordered = T,
                                    levels = c("scaffold1", "scaffold2", "scaffold3", "scaffold4", "scaffold5",
                                               "scaffold6", "scaffold7", "scaffold8", "scaffold9", "scaffold10"))
genome_selection_sites$pc <- factor(genome_selection_sites$pc , ordered = T,
                                    levels = c("PC1", "PC2", "PC3", "PC4", "PC5",
                                               "PC6", "PC7", "PC8", "PC9", "PC10"))

## Plot
genome_selection_plot <- genome_selection_sites %>%
  ggplot(aes(x=pos/10^6, y=neg_log_p_value)) +
  geom_point(size=0.2, alpha=0.5) +
  geom_smooth(color="purple", se =F) +
  geom_hline(yintercept = -log(0.05/dim(genome_selection)[1]), linetype = "dashed") + # This line marks a p-value cut off of 0.05 after Bonferroni correction
  theme_cowplot() +
  scale_x_continuous(breaks=seq(0, 100, 10)) +
  coord_cartesian(ylim=c(0, 40), expand = F) +
  xlab("position (Mbp)") +
  ylab("-log(p)") +
  facet_grid(pc~lg, scales="free_x", space="free_x") +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.title.x=element_text(),
        legend.position="none",
        text = element_text(size=10),
        axis.text = element_text(size=6))+
  theme(legend.key = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(legend.text=element_text(size=11)) +
  theme(panel.background = element_rect(fill = '#FAFAFA')) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
  theme(panel.border = element_blank())
ggsave(filename  = "/home/projects/dp_00007/people/hmon/Cult.vs.Wild/Cult.vs.Wild/GlobalVC_pcangsdSelection_10pcs.png", 
       plot=genome_selection_plot, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 300)
dev.off()



maf <- npyLoad(paste0(basedir, "/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30.maf.npy"))  
maf_df <- as.data.frame(maf)
write_tsv(maf_df, paste0(basedir, "/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde30.maf.tsv"), col_names = F)

#Get the annotation file 
cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.2.labels | awk '{split($0,a,"_"); print $1"\t"a[1]}' > /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.annot



#NGSadmix
# ngsTools +ngsadmix
module load perl/5.30.2         
module load samtools/1.11
module load ngstools/20190624
module load ngsadmix/32


for k in $(seq 1 5); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.beagle.gz -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 1000 -o $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_1kiterNGSadmix.$k
done
for k in $(seq 6 12); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.beagle.gz -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 1000 -o $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_1kiterNGSadmix.$k
done
