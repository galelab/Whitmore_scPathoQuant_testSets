#!/bin/sh


#Author: Leanne Whitmore 
# Downloading data for data which was shown to have a coinfection with SARS-CoV-2 and hMPV 
# Data from 2 papers, one that detected a co infection in human bronchial data and the paper to the data it did that with 
#-- 1) Viral track paper https://www.sciencedirect.com/science/article/pii/S0092867420305687
#-- 2) Data paper https://doi.org/10.1101/2020.02.23.20026690
## options 
# -X max size to download (units are in kb)
mkdir sra_data
cd sra_data


#Hevere covid patient 1
prefetch SRR11181956.1 -X 100000000
./fastq-dump SRR11181956.1 --split-files
wait

#Hevere covid patient 2
prefetch SRR11181958.1 -X 100000000
./fastq-dump SRR11181958.1 --split-files

#Hevere covid patient 3
prefetch SRR11181959.1 -X 100000000
./fastq-dump SRR11181959.1 --split-files
wait

