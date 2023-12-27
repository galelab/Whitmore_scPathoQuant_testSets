# scPathoQuant testSets 
This repository contains all codes for testing software [```scPathoQuant```](https://doi.org/10.1101/2023.07.21.549987) which can be download from [```here```](https://github.com/galelab/scPathoQuant)


### HIV test set (Bradley_et_al_testSet/)
-----------------------------------------
Data to test scPathoQuant used from the following paper Bradley T, Ferrari G, Haynes BF, Margolis DM, Browne EP. Single-Cell Analysis of Quiescent HIV Infection Reveals Host Transcriptional Profiles that Regulate Proviral Latency. Cell Rep. 2018 Oct 2;25(1):107-117.e3. doi: 10.1016/j.celrep.2018.09.020. PMID: 30282021; PMCID: PMC6258175.

Codes are run in the order of which they are numbered 
* 1.DownloadData.sh - Download raw data from SRA website
* 2a.SCprocessinghg38_HIV_sb.sh - align raw data to human genome with HIV genome integrated 
* 2b.SCprocessinghg38_sb.sh - align raw to human genome only 
* 3.scPathoQuant.sh - align results from 2b to HIV genome
* 4.SC_Analysis.r - run analysis on data generated from scPathoQuant (results for scPathoQuant paper)
* genomefiles/ - contains HIV genome files 

### Sars-Cov2 test set (Ravindra_et_al_testSet/)
------------------------------------------------
Data to test scPathoQuant used from the following paper Ravindra NG, Alfajaro MM, Gasque V, Huston NC, Wan H, Szigeti-Buck K, Yasumoto Y, Greaney AM, Habet V, Chow RD, Chen JS, Wei J, Filler RB, Wang B, Wang G, Niklason LE, Montgomery RR, Eisenbarth SC, Chen S, Williams A, Iwasaki A, Horvath TL, Foxman EF, Pierce RW, Pyle AM, van Dijk D, Wilen CB. Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. PLoS Biol. 2021 Mar 17;19(3):e3001143. doi: 10.1371/journal.pbio.3001143. PMID: 33730024; PMCID: PMC8007021.

Codes are run in the order of which they are numbered 
* 1.DownloadData.sh - Download raw data from SRA website
* 2a.CellRanger_Sars-CoV2_sb.sh - align raw data to human genome with HIV genome integrated 
* 2b.CellRanger_sb.sh - align raw to human genome only 
* 3.scPathoQuant.sh - align results from 2b to HIV genome
* 4.SC_Analysis.r - run analysis on data generated from scPathoQuant (results for scPathoQuant paper)
* genomefiles/ - contains sars-cov2 genome files

 

 

