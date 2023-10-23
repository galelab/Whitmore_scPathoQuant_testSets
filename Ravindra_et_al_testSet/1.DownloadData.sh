#!/bin/sh

#Author: Leanne Whitmore 
# Downloading single cell data for cultured human bronchial epithelial cells infected with SARS-COV-19
# Data from paper 
# --Single-cell longitudinal analysis of SARS-CoV-2 infection in the human airay epithelium identifies target cells alterations in gene expression and cell state changes
# --SRA data https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA701930&o=acc_s%3Aa

#Note: I renamed folders after downloading manual so that SRR identifiers would be replaced with sample names

## options 
# -X max size to download (units are in kb)
mkdir sra_raw_data
cd sra_raw_data

#Control 
prefetch SRR13711613.1 -X 100000000
# /share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR13711613.1 --split-files
wait

#Co-cultured with SARS-CoV-2 for 1 day
prefetch SRR13711614.1 -X 100000000
# /share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR13711614.1 --split-files
wait

#Co-cultured with SARS-CoV-2 for 2 days
prefetch SRR13711615.1 -X 100000000
# /share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR13711615.1 --split-files
wait

#Co-cultured with SARS-CoV-2 for 3 days
prefetch SRR13711616.1 -X 100000000
# /share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR13711616.1 --split-files
wait