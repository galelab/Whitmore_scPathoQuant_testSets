#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -t 5800-12
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out



#AUTHOR: Leanne Whitmore
#DESCRIPTION: human bronch cells, cellranger v6


GENOME_DIR="/share/tools/10x/refdata-gex-GRCh38-2020-A/"
PATH_FASTQ="../sra_raw_data/"
samples=("$PATH_FASTQ"*)

# folder for analysis using scviralquant 
mkdir 2.CellRanger_GRCh38
cd 2.CellRanger_GRCh38

for folder in ${samples[@]}
do

    sample_name=${folder#$PATH_FASTQ}
    echo "Processing sample $sample_name"
    echo "$folder"

    srun -c 1 cellranger count --sample=$sample_name --id="$sample_name" \
    --jobmode=slurm --transcriptome=$GENOME_DIR --fastqs=$folder &
    wait

done

