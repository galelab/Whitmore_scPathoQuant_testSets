#!/bin/sh
#Author Leanne  Whitmore 
#Description run scViralQuant (download from here https://github.com/galelab/scViralQuant)

PATH_FASTQ="./2.CellRanger_hg38/"
samples=("$PATH_FASTQ"SRR6825*)

#bowtie2

mkdir scViralQuantResuls_bowtie2
for folder in ${samples[@]}
do
    echo "$folder"
    sample_name=${folder#$PATH_FASTQ}
    mkdir scViralQuantResuls_bowtie2/"$sample_name"
    echo "Processing sample $sample_name"
    scviralquant.py -10x "$folder" --alignment_type="global" --aligner="bowtie2" -p 16 \
    -p2genome genomefiles -op scViralQuantResuls_bowtie2/"$sample_name"
    wait

done

mkdir scViralQuantResuls_bbmap
for folder in ${samples[@]}
do
    sample_name=${folder#$PATH_FASTQ}
    mkdir scViralQuantResuls_bbmap/"$sample_name"
    echo "Processing sample $sample_name"
    scviralquant.py -10x "$folder" --alignment_type="global" --aligner="bbmap" -p 16 \
    -p2genome genomefiles -op scViralQuantResuls_bbmap/"$sample_name"
    wait

done
