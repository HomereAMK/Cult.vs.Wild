#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/Cult.vs.Wild
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N LDApr23
#PBS -e LDApr23.err
#PBS -o LDApr23.out
#PBS -l nodes=2:ppn=20:thinnode
#PBS -l walltime=200:00:00
#PBS -l mem=100gb
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

#Load module for R
module load gsl/2.6
module load imagemagick/7.0.10-13
module load gdal/2.2.3
module load geos/3.8.0
module load jags/4.2.0
module load hdf5
module load netcdf
module load boost/1.74.0
module load openssl/1.0.0
module load lapack
module load udunits/2.2.26
module load proj/7.0.0
module load gcc/10.2.0
module load intel/perflibs/64/2020_update2
module load R/4.0.0



# LD pruning
##### Gets .pos file:
TRIAL=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221
zcat $TRIAL.mafs.gz | tail -n +2 | cut -f 1,2 > $TRIAL.LD.pos
##### Run ngsLD (needs to give the number of sites)
NSITES=`zcat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.pos.gz | wc -l `
/home/projects/dp_00007/apps/ngsLD/ngsLD --n_threads 40 --geno $TRIAL.beagle.gz --probs --n_ind 100 --n_sites $NSITES --pos $TRIAL.LD.pos --max_kb_dist 20 | pigz -p 40 >  $TRIAL.max_kb_dist20.LD.gz
##### Set variables-- the chromosomes
CHRs=("scaffold1" "scaffold2" "scaffold3" "scaffold4" "scaffold5" "scaffold6" "scaffold7" "scaffold8" "scaffold9" "scaffold10")
##### Splits ngsLD file per chromosome:
for query in ${CHRs[*]}
    do
    zcat $TRIAL.max_kb_dist20.LD.gz | grep "${query}" | pigz -p 40 > $TRIAL.max_kb_dist20.LD."${query}".gz
done
wait
##### Get the LD files list
find $TRIAL.max_kb_dist20.LD.*.gz > $TRIAL.max_kb_dist20.LD.PerCHR.list
##### LD pruning:
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput
parallel --will-cite --dryrun python /home/projects/dp_00007/people/hmon/EUostrea/00_scripts/prune_ngsLD.py --input {} --max_dist 100000 --min_weight 0.2 --output  $OUTPUTFOLDER/{/.}.17apr23.pruneIN --print_excl  $OUTPUTFOLDER/{/.}.17jan23.pruneOUT :::: $TRIAL.max_kb_dist20.LD.PerCHR.list > $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.RunsLDpruning.txt
##### Create the .sh
for num in `seq 0 10`
do
    num1=`expr $num + 1`
    # echo $num1
    head -n $num1 $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.RunsLDpruning.txt | tail -n1 > $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.SplitLD_RunsLDpruningUp_${num}.sh
done
wait
##### chmod the .sh
for num in `seq 0 10`
do
    chmod +x $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.SplitLD_RunsLDpruningUp_${num}.sh
done
wait
#### 
for num in `seq 0 10`
do
$TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.SplitLD_RunsLDpruningUp_${num}.sh
done
wait

cat $OUTPUTFOLDER/*.pruneIN > $OUTPUTFOLDER/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.pruneIN

# rm(list=ls(all.names = TRUE))
# library(tidyverse)
# basedir="/home/projects/dp_00007/people/hmon/Cult.vs.Wild"  
# chr<- c("scaffold1:", "scaffold2:", "scaffold3:", "scaffold4:", "scaffold5:",
#         "scaffold6:", "scaffold7:", "scaffold8:", "scaffold9:", "scaffold10:")
# pruned_position <- read_lines(paste0(basedir, "/02_angsdOutput/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.pruneIN")) %>% 
#   str_remove("scaffold10:") %>% str_remove("scaffold1:") %>% str_remove("scaffold2:") %>% 
#   str_remove("scaffold3:") %>% str_remove("scaffold4:") %>% str_remove("scaffold5:") %>%
#   str_remove("scaffold6:") %>% str_remove("scaffold7:") %>% str_remove("scaffold8:") %>%
#   str_remove("scaffold9:") %>% 
#   as.integer()
# pruned_position
# pruned_snp_list <- read_tsv(paste0(basedir, "/02_angsdOutput/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.mafs.gz"))%>%
#   dplyr::select(1:4) %>%
#   filter(position %in% pruned_position)
# write_tsv(pruned_snp_list, paste0(basedir, "/02_angsdOutput/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt"), col_names = F)

angsd sites index $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt