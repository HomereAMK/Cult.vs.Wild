Aquaculture vs Wild: two comparative case study -----------
Loch Nell vs Loch Ryan and Sydkoster vs Svallhagen_Tjärnö(or RAMS, LILL, KALV)


> A new dataset which 
>
## Launcher
```bash
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

```
> Fst-based analysis
## get saf per pop with minIND= 25%
```bash
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
```bash
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
```bash
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
```bash
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

```bash
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
```bash
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
```

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

> I need the read depth before doing the variant calling
```bash
BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput
N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist | wc -l`
angsd -bam $BAMLIST \
-ref $REF \
-out $OUTPUTFOLDER/Apr23_baq_mincov_minq20_minmapq20 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-doDepth 1 -doCounts 1 -maxDepth 2000 -dumpCounts 1 \
-P 8 \
-minMapQ 20 -minQ 20 \
-remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 
```

> On c2 get the .presenceGlobal
```R
nind_angsd <- data.table::fread("02_angsdOutput/qsubApr23_baq_mincov_minq20_minmapq20.mafs.gz")%>%
  count(nInd)
write_tsv(nind_angsd, "02_angsdOutput/qsubApr23_baq_mincov_minq20_minmapq20.presenceGlobal")
```

> Angsd Variant calling
```bash
BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput
N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist | wc -l`
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'
MINDP=73 # Minimum depth filter
MAXDP=221 # Maximum depth filter
MINQ=20 # Minimum quality filter
MINMAF=0.01 # Minimum minor allele frequency filter
MINMAPQ=20 # Minimum mapping quality (alignment score) filter, default value is 20

  angsd \
  -bam $BAMLIST \
  -ref $REF \
  -out $OUTPUTFOLDER/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221 \
  -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
  -setMinDepth $MINDP -setMaxDepth $MAXDP \
  -minQ $MINQ -minMapQ $MINMAPQ \
  -SNP_pval 1e-6 -minMaf $MINMAF \
  -nThreads 40 \
  $EXTRA_ARG

## Create a SNP list to use in downstream analyses
gunzip -c $OUTPUTFOLDER/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.mafs.gz | cut -f 1,2,3,4 | tail -n +2 > $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt

angsd sites index $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt

## Also make it in regions format for downstream analyses
cut -f 1,2 $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt | sed 's/\t/:/g' > $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_regions.txt

## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
cut -f1 $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt | sort | uniq > $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.chrs

```

>LD pruning for Admixture and PCA
##### Gets .pos file:
```
TRIAL=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221
zcat $TRIAL.mafs.gz | tail -n +2 | cut -f 1,2 > $TRIAL.LD.pos
```
##### Run ngsLD (needs to give the number of sites)
```bash
NSITES=`zcat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/02_angsdOutput/Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.pos.gz | wc -l `
/home/projects/dp_00007/apps/ngsLD/ngsLD --n_threads 40 --geno $TRIAL.beagle.gz --probs --n_ind 100 --n_sites $NSITES --pos $TRIAL.LD.pos --max_kb_dist 20 | pigz -p 40 >  $TRIAL.max_kb_dist20.LD.gz
```
##### Set variables-- the chromosomes
```
CHRs=("scaffold1" "scaffold2" "scaffold3" "scaffold4" "scaffold5" "scaffold6" "scaffold7" "scaffold8" "scaffold9" "scaffold10")
```

##### Splits ngsLD file per chromosome:
```
for query in ${CHRs[*]}
    do
    zcat $TRIAL.max_kb_dist20.LD.gz | grep "${query}" | pigz -p 40 > $TRIAL.max_kb_dist20.LD."${query}".gz
done
```
🤝
##### Get the LD files list
```
find $TRIAL.max_kb_dist20.LD.*.gz > $TRIAL.max_kb_dist20.LD.PerCHR.list
```
##### LD pruning:
```bash
parallel --will-cite --dryrun python $BASEDIR/00_scripts/prune_ngsLD.py --input {} --max_dist 100000 --min_weight 0.2 --output  $OUTPUTFOLDER/{/.}.17apr23.pruneIN --print_excl  $OUTPUTFOLDER/{/.}.17jan23.pruneOUT :::: $TRIAL.max_kb_dist20.LD.PerCHR.list > $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.RunsLDpruning.txt
```

🤝
```bash
for num in `seq 0 10`
do
    num1=`expr $num + 1`
    # echo $num1
    head -n $num1 $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.RunsLDpruning.txt | tail -n1 > $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.SplitLD_RunsLDpruningUp_${num}.sh
done
```


```bash
for num in `seq 0 10`
do
    chmod +x $TRIAL.max_kb_dist20kb.-min_weight0.2.17apr23.SplitLD_RunsLDpruningUp_${num}.sh
done
```

```bash
#will only launch the next job if the precedent is finished
job1.sh && job2.sh && job3.sh
python prune_ngsLD.py --input /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/VariantCalling/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.max_kb_dist100.LD.scaffold1.gz --max_dist 100000 --min_weight 0.2 --output /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.max_kb_dist100.LD.scaffold1.pruneIN --print_excl /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.max_kb_dist100.LD.scaffold1.pruneOUT
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
POP=("NELL" "RYAN" "TRAL" "SYDK" "RAMS" "KALV")
SITES=
for i1 in `seq 0 $((${#POP[@]}-2))`
do
    N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list | wc -l`
    angsd -b /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list \
    -anc $REF \
    -ref $REF \
    -out /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/17apr23--mindInd0.25_Unfolded_Cult.vs.Wild_${POP[i1]} \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 \
    -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) \
    -GL 1 -doSaf 1 -doGlf 3 -doMaf 1 -minmaf 0.05 -doHWE 1 -doMajorMinor 3 \
    -sites $SITES \
    -dosnpstat 1 -hwe_pval 1e-6 -SNP_pval 1e-6 
done


# mafs files are in /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/
#Jan23--Unfolded_Cult.vs.Wild_NELL.mafs.gz 
OUTPUTFOLDER2=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/04_Relatedness
POP=("NELL" "RYAN" "TRAL" "SYDK" "RAMS" "KALV")
for query in ${POP[*]}
do
  zcat /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/17apr23--mindInd0.25_Unfolded_Cult.vs.Wild_${query}.mafs.gz | cut -f6 | sed 1d > $OUTPUTFOLDER2/17apr23--mindInd0.25_Unfolded_Cult.vs.Wild_${query}.freq
done


for query in ${POP[*]}
do
N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list | wc -l`
$NGSRELATE \
-g /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/17apr23--mindInd0.25_Unfolded_Cult.vs.Wild_${query}.glf.gz \
-f $OUTPUTFOLDER2/17apr23--mindInd0.25_Unfolded_Cult.vs.Wild_${query}.freq \
-n $((N_IND*1)) -O $OUTPUTFOLDER2/17apr23--mindInd0.25_Unfolded_Cult.vs.Wild_${query}.res \
-p 10
done

```




> 19 apr23, Iqsub the angsd and pcangsd with the LDpruned list
```bash
BAMLIST=/home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Apr23--Cult.vs.Wild_VC_bamlist
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
THREADS=8
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
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub \
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
gunzip -c $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub.mafs.gz | cut -f 1,2,3,4 | tail -n +2 > $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub.txt

angsd sites index $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub.txt

## Also make it in regions format for downstream analyses
cut -f 1,2 $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub.txt | sed 's/\t/:/g' > $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub_regions.txt

## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
cut -f1 $OUTPUTFOLDER/global_snp_list_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221.txt | sort | uniq > $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub.chrs

wait

#pcangsd
pcangsd -b $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_iqsub.beagle.gz \
--selection \
--selection_e 10 \
--pcadapt \
--sites_save \
--snp_weights \
--pcadapt \
--maf_save \
--threads 30 \
--admix \
-e 10 \
-o $OUTPUTFOLDER/LDprunedlist_Apr23_VC_minq20_minMaf0.01_nominInd_setMinDepth73_setMaxDepth221_pcangsde10
```