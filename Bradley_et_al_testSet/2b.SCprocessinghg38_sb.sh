#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

#AUTHOR: Leanne Whitmore
#DESCRIPTION: Duke HIV cells, cellranger v6.1.2

mkdir /2.CellRanger_hg38
cd ./2.CellRanger_hg38
GENOME_DIR="./refdata-cellranger-hg38-pEGFP-N1-4.0.0"

PATH_FASTQ="./sra_data/"
samples=("$PATH_FASTQ"SRR6825026*)
echo "Number of samples = ${#samples[*]} (should be 4)"

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
