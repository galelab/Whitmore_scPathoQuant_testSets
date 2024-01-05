PATH=$PATH:/share/tools/10x/cellranger-6.1.2
export PATH
GENOME_DIR="./refdata-gex-GRCh38-2020-A/"
PATH_FASTQ="./sra_data/"
samples=("$PATH_FASTQ"SRR1118195*)

echo "Number of samples = ${#samples[*]} (should be 3)"

mkdir coinfectionresults
cd coinfectionresults

mkdir SPQresults

for folder in ${samples[@]}
do

    sample_name=${folder#$PATH_FASTQ}
    echo "Processing sample $sample_name"
    echo "Expected cell count ${expectedcells[counter]}"
    echo "$folder"

    cellranger count --sample=$sample_name --id="$sample_name" --jobmode=local \
        --localcores 16 --transcriptome=$GENOME_DIR --fastqs=$folder &
    wait

    mkdir SPQresults/"$sample_name"
    scpathoquant -10x "./coinfectionresults/"$sample_name --aligner="bbmap" -p 16 \
        -p2genome ./genomefiles -op SPQresults/"$sample_name" 
    wait
done
