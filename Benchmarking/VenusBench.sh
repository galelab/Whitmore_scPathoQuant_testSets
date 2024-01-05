repo_dir=./Venus
out_dir=./VenusBenchTrackResults

conda create --name venus --file ${repo_dir}/reference_files/venus_spec-file.txt
conda activate venus

python3.9 ${repo_dir}/src/module-index/module-index.py \
    --humanFASTA ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusFASTA ./genomefiles/pNL4-3renamed_HIVonly.fa \
    --module detection \
    --thread 8 \
    --out ${out_dir}/indices \
    --hGenome ${out_dir}/indices/human2.genomeDir \
    --virusGenome ${out_dir}/indices/HIV2.genomeDir

python3.9 ${repo_dir}/src/module-detection/module-detection.py \
    --read ./sra_data/SRR6825025/SRR6825025_S1_L001_R2_001.fastq ./sra_data/SRR6825025/SRR6825025_S1_L001_R1_001.fastq \
    --virusThreshold 1 \
    --virusChrRef ${repo_dir}/reference_files/virus_chr-ref_test.tsv \
    --virusGenome ${out_dir}/indices/HIV2.genomeDir \
    --humanGenome ${out_dir}/indices/human2.genomeDir \
    --readFilesCommand zcat \
    --thread 8 \
    --singleCellBarcode 1 16 \
    --singleUniqueMolIdent 17 10 \
    --singleWhitelist ${repo_dir}/reference_files/3M-february-2018.txt \
    --out ${out_dir}/detection/single-cell