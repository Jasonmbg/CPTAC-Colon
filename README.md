# CPTAC-Colon
Part II of the Athens Comprehensive Cancer Center Project (https://accc.gr/index.html) 

#### NOTE: *Core* parts of this REPOSITORY are also part of the broader ACCC-CRC project related to "*Molecular and functional profiling unravels targetable vulnerabilities in colorectal cancer*" accompanying manuscript.

#### Efstathios-Iason Vlachavas
###### DKFZ-Division of Molecular Genome Analysis (B050)
###### Efstathios-Iason.Vlachavas@dkfz-heidelberg.de
###### svlachavas@eie.gr

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10959700.svg)](https://doi.org/10.5281/zenodo.10959700)

Relative public dataset: coad_cptac_2019

Publication: Vasaikar et al. 2019 Cell. {https://www.ncbi.nlm.nih.gov/pubmed/31031003}

Data were downloaded from the open-source resource cBioPortal web portal (https://www.cbioportal.org/)
and the respective files were utilized (https://www.cbioportal.org/study/summary?id=coad_cptac_2019):

- File 1: **data_RNA_Seq_v2_expression_median.txt** which then saved as an .rda file for downstream analysis.

- File 2: **data_phosphoprotein_quantification.txt** which then saved as an .rda file for downstream analysis.

- File 3: **data_mutations_extended.txt** which then was further processed, filtered and saved as maf object and described in the respective Zenodo folder materials.

- File 4: **data_clinical_sample.txt**. Here we modified it to skip the first 4 # header lines and saved it initially as **coad_cptac_2019_clinical_data.tsv**. Then, after further processing and sample selection for analysis purposes, we ended we 2 phenotype-like objects, namely **CPTAC_clinical_start_INFO.tsv** and **CPTAC_clinical_final_mod_INFO.tsv**. These are the files that are used in the analyses reproduced in the .qmd files.

###### These initial files, along with modified-processed-additional results, are all included in the *Data* & *Mut_Freq* folders

##### **Note**: Any 'save' function in R (as well as other commonly used functions in R to export objects in various formats, such as *write_tsv*) documented along with a hashtag # sign are typically used to demonstrate how an object was created and is not intended to be run as part of regular code execution.

## Description

This public dataset was used as an additional proxy to increase sample size regarding the analysis of the *in-house* ACCC-CRC cohort, as well as to utilize the additional phosphoproteomic layer to provide further molecular insights. On this premise, the pathway relative activities that can be also inferred from the phosphoproteome data, can support putative mechanisms and molecular players extracted from the gene expression data.

## Notes on data acquisition and cohort definition

Metadata utilization
- Clinical data was directly utilized from the cBioPortal web page: further patient selection was performed as described in the manuscript, based on the presence of specific mutations. The relative files that will be used are **CPTAC_clinical_start_INFO.tsv** & mainly for downstream purposes **CPTAC_clinical_final_mod_INFO.tsv**.

RNASeq data
- The analysis starts essentially with the *semi-processed;RSEM estimated counts* 
**(RNA Seq V2 RSEM UQ Log2)**, gene annotation and sample phenotype information.

Phosproproteomics data
- Likewise, the bioinformatics workflow begins with the processed phosphoproteomics intensities (normalized log2 abundance ratios)

- Further methodological details on the complete analysis are included in the accompanied.qmd files & Materials & Methods section of the manuscript.

###### For file size, reproducibility and memory considerations, for both datasets we start with .rda files of the raw data versions.

Further details on each of the utilized databases/tools/repositories based on their respective publications:

- Omnipath: https://doi.org/10.15252/msb.20209923

- PTMSigDB: https://doi.org/10.1074/mcp.TIR118.000943

- decoupleR: https://doi.org/10.1093/bioadv/vbac016

- cBioPortal: https://doi.org/10.1158/2159-8290.CD-12-0095

## Neccesary R packages that need to be installed for reproducing the analysis for the CPTAC dataset:

```r

packages = c(
    "edgeR",
    "circlize",
    "ComplexHeatmap"
    "tidyverse",
    "data.table",
    "progeny",
    "xml2",
    "downlit",
    "here",
    "decoupleR",
    "OmnipathR",
    "limma",
    "PhosR",
    "AnnotationDbi"
    "org.Hs.eg.db",
    "ggplot2",
    "forcats",
    "DOSE",
    "clusterProfiler",
    "ggplot2"
)

if(!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
}

library(BiocManager)

for(p in packages) {
    if(!requireNamespace(p)) {
        install(p)
    }
}

```
## Implementation

- The user just needs to download/clone the respective github repository;
- For example `git clone https://github.com/Jasonmbg/CPTAC-Colon.git
