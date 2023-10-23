#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -t 5800-12
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

#AUTHOR: Leanne Whitmore
#DESCRIPTION: human bronch cells, cellranger v6

mkdir 2.CellRanger_GRCh38_SARS_CoV
cd 2.CellRanger_GRCh38_SARS_CoV
counter=0
GENOME_DIR="./refdata-cellranger-GRCh38_SarsCoV2/"

for folder in ${samples[@]}
do


    sample_name=${folder#$PATH_FASTQ}
    echo "Processing sample $sample_name"
    echo "Expected cell count ${expectedcells[counter]}"
    echo "$folder"

    srun -c 1 cellranger count --sample=$sample_name --id="$sample_name" --jobmode=slurm \
    --transcriptome=$GENOME_DIR --fastqs=$folder &

    wait
    counter=$((counter+1))

done
