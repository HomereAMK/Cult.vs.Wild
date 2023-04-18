#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/Cult.vs.Wild
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N pcangsdApr23
#PBS -e pcangsdApr23.err
#PBS -o pcangsdApr23.out
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

BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
THREADS=15
SNP_LIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'



## Re-run angsd LD Pruned SNPs list minweight0.2 Global dataset 
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221 \
-doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-minQ 20 -minMapQ 20 \
-P $THREADS \
$EXTRA_ARG \
-sites $SNP_LIST \
-rf $LG_LIST

#pcangsd
pcangsd -b $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.beagle.gz \
--selection \
--selection_e 10 \
--pcadapt \
--sites_save \
--snp_weights \
--pcadapt \
--maf_save \
--threads 40 \
--admix \
-e 10 \
-o $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde10
