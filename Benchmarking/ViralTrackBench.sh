
mkdir ViralTrackBenchResults 
cd ViralTrackBenchResults
# ###UMI tools 
Step 2: Identify correct cell barcodes
umi_tools whitelist --stdin ./sra_data/SRR6825025/SRR6825025_S1_L001_R1_001.fastq \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --set-cell-number=772 \
                    --log2stderr > whitelist.txt;
 # same number as 10x 
# Step 3: Extract barcdoes and UMIs and add to read names
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin ./sra_data/SRR6825025/SRR6825025_S1_L001_R1_001.fastq \
                  --stdout SRR6825025_S1_L001_R1_001.extracted.fastq.gz \
                  --read2-in ./sra_data/SRR6825025/SRR6825025_S1_L001_R2_001.fastq \
                  --read2-out=SRR6825025_S1_L001_R2_001.extracted.fastq.gz \
                  --whitelist=whitelist.txt; 

mkdir ../human_ref_genome/viraltrack_star_index
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ../human_ref_genome/viraltrack_star_index \
    --genomeFastaFiles ./genomefiles/pNL4-3renamed_HIVonly.fa ../human_ref_genome/genome.fa
wait 

Rscript ../tools/Viral-Track/Viral_Track_scanning.R ./Viral-Trackfiles/ParametersHIV.txt ./Viral-Trackfiles/targetfiles.txt
wait 
# mkdir ./HIVoutput/SRR6825025/
Rscript ./Viral-Track/Viral_Track_cell_demultiplexing.R ./Viral-Trackfiles/ParametersHIV.txt ./Viral-Trackfiles/targetfiles.txt
wait 