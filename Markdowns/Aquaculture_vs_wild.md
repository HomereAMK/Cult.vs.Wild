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

    POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV" "TRAL")

    for query in ${POP[*]}
    do 
        grep ${query} /home/projects/dp_00007/people/hmon/Flat_oysters/01_infofiles/bam_list_aug22.txt > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${query}-Fst.list

    done



## get maf/saf per pop with minIND= 10%
```
    REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
    POP=(NELL RYAN CLEW SYDK RAMS KALV TRAL)

for i1 in `seq 0 $((${#POP[@]}-2))`
do
    angsd -nThreads 10 \
	-b /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Jan23--Cult.vs.Wild_${POP[i1]}-Fst.list \
	-anc $REF \
	-ref $REF \
	-out /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Jan23--setMinDepth100Unfolded_Cult.vs.Wild_${POP[i1]} \
	-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 \
	-minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 100 -setMaxDepth 1200 \
	-remove_bads 1 -only_proper_pairs 1 -C 50 
done
```
🤝


## get the sfs step
    cd /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/
    POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV" "TRAL")
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

POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV" "TRAL")
for query in ${POP[*]}
    do 
    realSFS Jan23--setMinDepth100Unfolded_Cult.vs.Wild_NELL.saf.idx Jan23--setMinDepth100Unfolded_Cult.vs.Wild_${query}.saf.idx -fold 1 -P 40 > /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/NELL.${query}_setMinDepth100Jan23.sfs
done

## sfs with a sliding window step nolist 15kb 15kb
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






    realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/NELL.SYDK_Sept22.fst.idx -win 15000 -step 15000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-15000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst15kb/Sept22_Cult.vs.Wild_mindIND0.1_NELL.SYDK_15KB_15KB--Fst.tsv 



## sfs with a sliding window step nolist 5kb 5kb
#on the 4th oct 22 I run all the comparison except the ones that include KALV (saf file for KALV did not run for some reason)
    POP=("NELL" "RYAN" "CLEW" "SYDK" "RAMS" "KALV" "TRAL")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="${POP[i1]}"
            pop2="${POP[i2]}"
                realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_setMinDepth100Jan23.fst.idx -win 5000 -step 5000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-5000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst5kb/Jan23_Cult.vs.Wild_propersetMinDepth100Jan23_${POP[i1]}.${POP[i2]}_5KB_5KB--Fst.tsv 
        done
    done

