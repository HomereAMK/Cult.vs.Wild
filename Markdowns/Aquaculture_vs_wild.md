Aquaculture vs Wild: two comparative case study -----------
Loch Nell vs Loch Ryan and Sydkoster vs Svallhagen_Tjärnö(or RAMS, LILL, KALV)


> A new dataset which 
>
## Launcher
```
#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/Cult.vs.Wild
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N SFS.19jan23
#PBS -e SFS.19jan23.err
#PBS -o SFS.19jan23.out
#PBS -l nodes=1:ppn=10:thinnode
#PBS -l walltime=400:00:00
#PBS -l mem=80gb
#PBS -m n
#PBS -r n
```

## Angsd c2 module loading
    module load tools ngs computerome_utils/2.0
    module load htslib/1.16
    module load bedtools/2.30.0
    module load pigz/2.3.4
    module load parallel/20210722
    module load angsd/0.940



## Create fst list 

    POP=("NELL" "RYAN" "TRAL" "SYDK" "RAMS" "KALV")

    for query in ${POP[*]}
    do 
        grep ${query} /home/projects/dp_00007/people/hmon/Flat_oysters/01_infofiles/bam_list_aug22.txt > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list

    done



## get saf per pop with minIND= 25%
```
for i1 in `seq 0 $((${#POP[@]}-2))`
do
    N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list | wc -l`
    /home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 10 \ 
    -b /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list \
    -anc $REF \
    -ref $REF \
    -out /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_Cult.vs.Wild_${POP[i1]} \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 \
    -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) \
    -GL 1 -doSaf 1
done
```
🤝

## get the sfs step
```
    cd /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/
    POP=("NELL" "RYAN" "SYDK" "RAMS" "KALV" "TRAL")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="Jan23--setMinDepth100Unfolded_Cult.vs.Wild_${POP[i1]}"
            pop2="Jan23--setMinDepth100Unfolded_Cult.vs.Wild_${POP[i2]}"
            N_SITES=`realSFS print $pop1.saf.idx $pop2.saf.idx | wc -l`
            echo -ne "${POP[i1]}\t${POP[i2]}\t$N_SITES\t"
            if [[ $N_SITES == 0 ]]; then
                echo "NA"
            else
                realSFS $pop1.saf.idx $pop2.saf.idx -fold 1 -P 40 > /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_setMinDepth100Jan23.sfs
                realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_setMinDepth100Jan23.sfs -fold 1 -P 40 -fstout /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_setMinDepth100Jan23
                realSFS fst stats /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_setMinDepth100Jan23.fst.idx -P 40
            fi
        done
    done > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/Jan23_Cult.vs.Wild_proper_setMinDepth100_Fst.tsv
```
🤝
```
    cd /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/
    POP=("NELL" "RYAN" "SYDK" "RAMS" "KALV" "TRAL")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="26Jan23--mindInd0.25_Unfolded_Cult.vs.Wild_${POP[i1]}"
            pop2="26Jan23--mindInd0.25_Unfolded_Cult.vs.Wild_${POP[i2]}"
            N_SITES=`realSFS print $pop1.saf.idx $pop2.saf.idx | wc -l`
            echo -ne "${POP[i1]}\t${POP[i2]}\t$N_SITES\t"
            if [[ $N_SITES == 0 ]]; then
                echo "NA"
            else
                realSFS $pop1.saf.idx $pop2.saf.idx -fold 1 -P 40 > /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}.sfs
                realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}.sfs -fold 1 -P 40 -fstout /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}
                realSFS fst stats /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}.fst.idx -P 40
            fi
        done
    done > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/26Jan23--mindInd0.25_Unfolded_Fst.tsv
```

## sfs with a sliding window step nolist 15kb 15kb
```
#on the 4th oct 22 I run all the comparison except the ones that include KALV (saf file for KALV did not run for some reason)
    POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV" "TRAL")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="${POP[i1]}"
            pop2="${POP[i2]}"
                realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_setMinDepth100Jan23.fst.idx -win 15000 -step 15000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-15000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst15kb/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_${POP[i1]}.${POP[i2]}_15KB_15KB--Fst.tsv 
        done
    done
```
🤝

```
    POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV" "TRAL")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="${POP[i1]}"
            pop2="${POP[i2]}"
                realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}.fst.idx -win 15000 -step 15000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-15000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst15kb/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}_15KB_15KB--Fst.tsv 
        done
    done
```


## sfs with a sliding window step nolist 5kb 5kb
    POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV" "TRAL")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="${POP[i1]}"
            pop2="${POP[i2]}"
                realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}.fst.idx -win 5000 -step 5000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-5000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst5kb/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_${POP[i1]}.${POP[i2]}_5KB_5KB--Fst.tsv 
        done
    done
🤝

    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="${POP[i1]}"
            pop2="${POP[i2]}"
                realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}.fst.idx -win 5000 -step 5000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-5000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst5kb/26Jan23--mindInd0.25_Unfolded_${POP[i1]}.${POP[i2]}_5KB_5KB--Fst.tsv 
        done
    done

## sfs with a sliding window step nolist 1kb 1kb
    POP=("NELL" "RYAN" "SYDK" "RAMS" "KALV" "TRAL")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="${POP[i1]}"
            pop2="${POP[i2]}"
                realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_setMinDepth100Jan23.fst.idx -win 1000 -step 1000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-1000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst1kb/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_${POP[i1]}.${POP[i2]}_1KB_1KB--Fst.tsv 
        done
    done

🤝

## Mapped the fasta Fst outliers to the genome
```bash
REF2=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G
OTL=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/chr2_NELLvsRYANSYDK_fst_outlier.fasta

bowtie2-build /home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta /home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G
🤝

bowtie2 \
-r --very-sensitive \
-p 8 \
-x $REF2 \
-f $OTL
-F 1000,1000 \
-S /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/chr2_NELLvsRYANSYDK_fst_outlier_mapped_to_gadMor2.sam


bowtie2 -r --very-sensitive -p 8 -x $REF2 -f $OTL -S /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/chr2_NELLvsRYANSYDK_fst_outlier_mapped_to_gadMor2.sam

bowtie2 -r --very-sensitive -p 8 -x $REF2  -U $OTL -S /path/to/outputfile.sam


samtools view -q 20 /workdir/cod/gosl-cod/angsd/popminind5/bam_list_realigned_mincov_filtered_mindp60_maxdp337_minind33_minq20_popminind5_fst_outlier_mapped_to_gadMor2.sam\
> /workdir/cod/gosl-cod/angsd/popminind5/bam_list_realigned_mincov_filtered_mindp60_maxdp337_minind33_minq20_popminind5_fst_outlier_mapped_to_gadMor2_minmapq20.sam
```


## Variant calling
> Modules
```bash
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
```
> create a bamlist
```bash
POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV")
for query in ${POP[*]}
do 
    grep ${query} /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list >> /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
done
```
🤝
> Angsd Variant calling
```bash
BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput
N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist | wc -l`
  angsd \
  -bam $BAMLIST \
  -ref $REF \
  -out $OUTPUTFOLDER/Apr23_VC_minq20_minind25per_nomaxmindp \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
  -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) \
  -doCounts 1 -dumpCounts 2 \
  -GL 1 -doGlf 2 \
  -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 \
  -doIBS 1 -doCov 1 -makeMatrix 1 \
  -nThreads 40 
```
🤝
```
-> Output filenames:
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.arg"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.pos.gz"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.counts.gz"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.beagle.gz"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.mafs.gz"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.geno.gz"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.ibs.gz"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.ibsMat"
                ->"/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minind25per_nomaxmindp.covMat"
        -> Tue Apr 11 14:37:23 2023
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 832225206
        -> Number of sites retained after filtering: 8797487
        [ALL done] cpu-time used =  201445.69 sec
        [ALL done] walltime used =  71992.00 sec
```

>PCAngsd 
```bash
pcangsd -b $OUTPUTFOLDER/Apr23_VC_minq20_minind25per_nomaxmindp.beagle.gz \
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
-o $OUTPUTFOLDER/Apr23_VC_minq20_minind25per_nomaxmindp_pcangsde10
```

> Data wrangling on c2
```R
rm(list=ls(all.names = TRUE))
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)
basedir="/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput"  
genome_selection <- npyLoad(paste0(basedir, "/Apr23_VC_minq20_minind25per_nomaxmindp_pcangsde10.selection.npy"))  
genome_selection_df <- as.data.frame(genome_selection)
write_tsv(genome_selection_df, paste0(basedir, "/Apr23_VC_minq20_minind25per_nomaxmindp_pcangsde10.selection.tsv"), col_names = F)
```

> Relatedness
```bash
#get glf
BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV")
for i1 in `seq 0 $((${#POP[@]}-2))`
do
    N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list | wc -l`
    /home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 40 \ 
    -b /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list \
    -anc $REF \
    -ref $REF \
    -out /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/10apr23--mindInd0.25_Unfolded_Cult.vs.Wild_${POP[i1]} \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 \
    -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) \
    -GL 1 -doSaf 1 -doGlf 3 -doMaf 1 -minmaf 0.05 -doHWE 1 -doMajorMinor 3 \
    -dosnpstat 1 -hwe_pval 1e-6 -SNP_pval 1e-6 
done


# mafs files are in /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/
#Jan23--Unfolded_Cult.vs.Wild_NELL.mafs.gz 
OUTPUTFOLDER2=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/04_Relatedness
POP=("NELL" "RYAN" "TRAL" "SYDK" "RAMS" "KALV")
for query in ${POP[*]}
do
  zcat /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Jan23--Unfolded_Cult.vs.Wild_${query}.mafs.gz | cut -f6 | sed 1d > $OUTPUTFOLDER2/10apr23_ngsrelate_minind0.75_LDprunedList_${query}.freq
done


for query in ${POP[*]}
do
N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list | wc -l`
$NGSRELATE \
-g /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Jan23--Unfolded_Cult.vs.Wild_${query}.glf.gz \
-f $OUTPUTFOLDER2/10apr23_ngsrelate_minind0.75_LDprunedList_${query}.freq \
-n $((N_IND*1)) -O $OUTPUTFOLDER2/10apr23_ngsrelate_minind0.75_LDprunedList_${query}.res \
-p 10
done

```

