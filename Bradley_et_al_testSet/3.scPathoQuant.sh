#!/bin/sh
#Author Leanne  Whitmore 
#Description run scViralQuant (download from here https://github.com/galelab/scViralQuant)

PATH_FASTQ="./2.CellRanger_hg38/"
samples=("$PATH_FASTQ"SRR6825*)

#bowtie2

mkdir scPathoQuantResults_bowtie2
for folder in ${samples[@]}
do
    echo "$folder"
    sample_name=${folder#$PATH_FASTQ}
    mkdir scPathoQuantResults_bowtie2/"$sample_name"
    echo "Processing sample $sample_name"
    scpathoquant -10x "$folder" --aligner="bowtie2" -p 16 \
    -p2genome genomefiles -op scPathoQuantResults_bowtie2/"$sample_name"
    wait

done

mkdir scPathoQuantResults_bbmap
for folder in ${samples[@]}
do
    sample_name=${folder#$PATH_FASTQ}
    mkdir scPathoQuantResults_bbmap/"$sample_name"
    echo "Processing sample $sample_name"
    scpathoquant -10x "$folder" --aligner="bbmap" -p 16 \
    -p2genome genomefiles -op scPathoQuantResults_bbmap/"$sample_name"
    wait

done
