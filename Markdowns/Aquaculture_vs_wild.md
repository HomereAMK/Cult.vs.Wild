Aquaculture vs Wild: two comparative case study -----------
Loch Nell vs Loch Ryan and Sydkoster vs Svallhagen_Tjärnö(or RAMS, LILL, KALV)


> A new dataset which 
## Angsd c2 module loading
    module load tools ngs computerome_utils/2.0
    module load htslib/1.16
    module load bedtools/2.30.0
    module load pigz/2.3.4
    module load parallel/20210722
    module load angsd/0.937



## Create fst list 

    POP=("NELL" "RYAN" "SYDK" "SVAL" "RAMS" "LILL" "KALV")

    for query in ${POP[*]}
    do 
        grep ${query} /home/projects/dp_00007/people/hmon/Flat_oysters/01_infofiles/bam_list_aug22.txt > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Sept22--Cult.vs.Wild_${query}-Fst.list

    done



## get maf per pop with minIND= 10%
    REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
    POP=(NELL RYAN SYDK SVAL RAMS LILL KALV)
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        pop1="${POP[i1]}"
        N_IND=`cat /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Sept22--Cult.vs.Wild_${POP[i1]}-Fst.list | wc -l`
        /home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 40 -ref $REF -anc $REF -bam /home/projects/dp_00007/people/hmon/Cult.vs.Wild/01_infofiles/Sept22--Cult.vs.Wild_${POP[i1]}-Fst.list -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 20 -minQ 20 -minInd $((N_IND1/10)) -GL 1 -doSaf 1 -out /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/Sept22--Unfolded_Cult.vs.Wild_${POP[i1]}
        done
    done
$$$$$$$$$$$$$$$$$$
# get the sfs step
    cd /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild
    POP=("NELL" "RYAN" "SYDK" "SVAL" "RAMS" "LILL" "KALV")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="Sept22--Unfolded_Cult.vs.Wild_${POP[i1]}"
            pop2="Sept22--Unfolded_Cult.vs.Wild_${POP[i2]}"
            N_SITES=`realSFS print $pop1.saf.idx $pop2.saf.idx | wc -l`
            echo -ne "${POP[i1]}\t${POP[i2]}\t$N_SITES\t"
            if [[ $N_SITES == 0 ]]; then
                echo "NA"
            else
                realSFS $pop1.saf.idx $pop2.saf.idx -fold 1 -P 40 > /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_Sept22.sfs
                realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_Sept22.sfs -fold 1 -P 40 -fstout /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_Sept22
                realSFS fst stats /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_Sept22.fst.idx -P 40
            fi
        done
    done > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/Sept22_Cult.vs.Wild_mindIND0.1--Fst.tsv


## sfs with a sliding window step nolist 15kb 15kb
    POP=("NELL" "RYAN" "SYDK" "SVALL" "RAMS" "LILL" "KALV")
    for i1 in `seq 0 $((${#POP[@]}-2))`
    do
        for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
        do
            pop1="${POP[i1]}"
            pop2="${POP[i2]}"
                realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/Cult.vs.Wild/${POP[i1]}.${POP[i2]}_Sept22.fst.idx -win 15000 -step 15000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-15000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/people/hmon/Cult.vs.Wild/03_Fst_based/SLFst15kb/Sept22_Cult.vs.Wild_mindIND0.1_${POP[i1]}.${POP[i2]}_15KB_15KB--Fst.tsv 
        done
    done