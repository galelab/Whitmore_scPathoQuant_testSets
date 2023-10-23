#!/bin/sh
# Author: Leanne Whitmore 
# Downloading data for CD4 HIV infected cells uses GFP as a marker of cell infection. Will use this as a positive 
# control to test scViralQuant
# Data from paper 
# --Single-cell RNA-seq of a primary CD4 T cell model of human HIV latency
# Bradley T, Ferrari G, Haynes BF, Margolis DM, Browne EP. Single-Cell Analysis of Quiescent HIV Infection Reveals Host Transcriptional Profiles that Regulate Proviral Latency. Cell Rep. 2018 Oct 2;25(1):107-117.e3. doi: 10.1016/j.celrep.2018.09.020. PMID: 30282021; PMCID: PMC6258175.

## options 
# -X max size to download (units are in kb)
mkdir sra_data
cd sra_data

#HIV INFECTED REP 1
# prefetch SRR6825024.1 -X 100000000
/share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR6825024.1 --split-files
wait
#HIV INFECTED REP 2
# prefetch SRR6825025.1 -X 100000000
/share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR6825025.1 --split-files
wait

#HIV UNINFECTED REP 1
# prefetch SRR6825026.1 -X 100000000
/share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR6825026.1 --split-files
wait
#HIV UNINFECTED REP 2
# prefetch SRR6825027.1 -X 100000000
/share/tools/sratoolkit.2.10.4-centos_linux64/bin/fastq-dump SRR6825027.1 --split-files
wait