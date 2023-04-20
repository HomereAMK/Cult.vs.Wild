#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/Cult.vs.Wild
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N FstApr23
#PBS -e FstApr23.err
#PBS -o FstApr23.out
#PBS -l nodes=2:ppn=20:thinnode
#PBS -l walltime=200:00:00
#PBS -l mem=100gb
#PBS -m n
#PBS -r n

## Angsd c2 module loading
    module load tools ngs computerome_utils/2.0
    module load htslib/1.16
    module load bedtools/2.30.0
    module load pigz/2.3.4
    module load parallel/20210722
    module load angsd/0.940


BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
GLOBALSNPLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt
POP=("NELL" "RYAN" "TRAL" "SYDK" "RAMS" "KALV")
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'
THREADS=10


for query in ${POP[*]}
do 
grep ${query} /home/projects/dp_00007/people/hmon/Flat_oysters/01_infofiles/bam_list_aug22.txt > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list
done
wait 

# get saf per pop with minIND= 25%
for i1 in `seq 0 $((${#POP[@]}-2))`
do
    N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list | wc -l`
    angsd -b /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list \
    -anc $REF \
    -ref $REF \
    -out /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]} \
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
    -P $THREADS \
    -minInd $((N_IND*1/4)) -minQ 20 -minMapQ 20 \
    -sites $GLOBALSNPLIST -rf $LG_LIST $EXTRA_ARG
done

wait 

# Estimate Fst in each pair of populations
cd  /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/
POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV")

for i1 in `seq 0 $((${#POP[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
    do
        pop1="Apr23_globallist_mindInd0.25_${POP[i1]}"
        pop2="Apr23_globallist_mindInd0.25_${POP[i2]}"
        N_SITES=`realSFS print $pop1.saf.idx $pop2.saf.idx | wc -l`
        realSFS $pop1.saf.idx $pop2.saf.idx -fold 1 -P 40 > /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}.sfs
        realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}.sfs -fold 1 -P 40 -fstout /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}
        realSFS fst stats /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}.fst.idx -P 40      
    done
done > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/Apr23_globallist_mindInd0.25_Fst_21apr23.tsv

wait

for i1 in `seq 0 $((${#POP[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
     do
        pop1="${POP[i1]}"
        pop2="${POP[i2]}"
        realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}.fst.idx -win 15000 -step 15000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-15000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}_15KB_15KB--Fst2.tsv 
    done
done

wait

for i1 in `seq 0 $((${#POP[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
     do
        pop1="${POP[i1]}"
        pop2="${POP[i2]}"
        realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}.fst.idx -win 1000 -step 1000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-1000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Apr23_globallist_mindInd0.25_${POP[i1]}.${POP[i2]}_1KB_1KB--Fst.tsv 
    done
done
