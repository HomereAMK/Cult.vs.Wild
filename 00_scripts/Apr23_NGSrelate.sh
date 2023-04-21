#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/Cult.vs.Wild
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N NGSrelateApr23
#PBS -e NGSrelateApr23.err
#PBS -o NGSrelateApr23.out
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


#get glf
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
OP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV")
SNP_LIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt
NGSRELATE=/home/projects/dp_00007/people/hmon/Software/ngsRelate/ngsRelate
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
THREADS=12

for query in ${POP[*]}
do
    N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list | wc -l`
    angsd -bam /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list \
    -ref $REF \
    -P $THREADS \
    -out /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_LDprunedSNPlist_mindInd0.75_Unfolded_Cult.vs.Wild_${query} \
    -GL 1 -doGlf 3 -doMaf 1 -minmaf 0.05 -doMajorMinor 3 \
    -sites $SNP_LIST -rf $LG_LIST \
    -minMapQ 20 -minQ 20 \
    -minInd $((N_IND*2/3)) -SNP_pval 1e-6 
done

wait

# mafs files are in /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/
#Jan23--Unfolded_Cult.vs.Wild_NELL.mafs.gz 
OUTPUTFOLDER2=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/04_Relatedness
POP=("NELL" "RYAN" "TRAL" "SYDK" "RAMS" "KALV")
for query in ${POP[*]}
do
  zcat /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_LDprunedSNPlist_mindInd0.75_Unfolded_Cult.vs.Wild_${query}.mafs.gz | cut -f6 | sed 1d > $OUTPUTFOLDER2/Apr23_LDprunedSNPlist_mindInd0.75_Unfolded_Cult.vs.Wild_${query}.freq
done

wait
for query in ${POP[*]}
do
N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list | wc -l`
$NGSRELATE \
-g /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_LDprunedSNPlist_mindInd0.75_Unfolded_Cult.vs.Wild_${query}.glf.gz \
-f $OUTPUTFOLDER2/Apr23_LDprunedSNPlist_mindInd0.75_Unfolded_Cult.vs.Wild_${query}.freq \
-n $((N_IND*1)) -O $OUTPUTFOLDER2/Apr23_LDprunedSNPlist_mindInd0.75_Unfolded_Cult.vs.Wild_${query}.res \
-p 10
done

