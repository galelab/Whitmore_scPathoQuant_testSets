# will use CD4_HIV_DukeData for benchmarking
mkdir PathogenTrackBenchResults 
cd PathogenTrackBenchResults

gzip -d barcodes.tsv.gz

mkdir ../human_ref_genome/pathogentrack_star_index
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ../human_ref_genome/pathogentrack_star_index \
     --genomeFastaFiles ../human_ref_genome/genome.fa \
     --sjdbGTFfile ../human_ref_genome/genes.gtf \
     --sjdbOverhang 100

python ./PathogenTrack/PathogenTrack.py count --project_id SRR6825025 \
        --star_index ../human_ref_genome/pathogentrack_star_index \
        --kraken_db ../tools/minikraken_8GB_20200312 \
        --barcode barcodes.tsv \
        --read1 ./sra_data/SRR6825025/SRR6825025_S1_L001_R1_001.fastq \
        --read2 ./sra_data/SRR6825025/SRR6825025_S1_L001_R2_001.fastq
