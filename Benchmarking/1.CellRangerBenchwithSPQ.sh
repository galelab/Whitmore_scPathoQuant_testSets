
PATH=$PATH:/share/tools/10x/cellranger-6.1.2
export PATH
GENOME_DIR="./refdata-cellranger-hg38-pEGFP-N1-4.0.0"
PATH_FASTQ="./sra_data/"
samples=("$PATH_FASTQ"SRR6825025*)

echo "Number of samples = ${#samples[*]} (should be 1)"
mkdir CellRangerwithscpathoquantBenchResults
cd CellRangerwithscpathoquantBenchResults

mkdir SPQResults
for folder in ${samples[@]}
do

    sample_name=${folder#$PATH_FASTQ}
    echo "Processing sample $sample_name"
    echo "Expected cell count ${expectedcells[counter]}"
    echo "$folder"

    cellranger count --sample=$sample_name --id="$sample_name" --jobmode=local \
        --localcores 8 --transcriptome=$GENOME_DIR --fastqs=$folder &

    wait
    mkdir SPQResults/"$sample_name"
    scpathoquant -10x ./CellRangerwithscpathoquantBenchResults/SRR6825025 \
        --aligner="bbmap" -p 8 \
        -p2genome ./genomefiles -op ./SPQResults/"$sample_name" 
    wait
    counter=$((counter+1))

done
