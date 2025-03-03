# RNAseq Snakemake Pipeline  

## Overview  

This repository contains a Snakemake pipeline for processing RNA-seq data. The pipeline automates quality control, read alignment, genes & isoforms quantification. It is designed to be scalable, reproducible, and easy to configure.  
## Features  

- **Automated Workflow:** Handles all major RNA-seq processing steps.  
- **Reproducibility:** Uses Snakemake for workflow management.  
- **Parallel Execution:** Supports multi-threading for faster processing.  
- **Modular Design:** Easily customizable for different organisms and experiments.  
- **Logging & Reports:** Generates QC and summary reports.  

## Workflow Steps  

1. **Quality Control** (`FastQC`)  
2. **Trimming & Filtering** (`Trim Galore`)  
3. **Read Alignment** (`STAR`)  
4. **Transcript Quantification** ( `Salmon`)
5. **Gene Quantification** (`featureCounts`)  
5. **Differential Expression Analysis** (`DESeq2`)    

## Installation  

### Prerequisites  

Ensure you have the following installed:  

- **Conda** (via Miniconda or Anaconda)  
- **Snakemake** (`conda install -c bioconda snakemake`)    
- **Required bioinformatics tools** (see `envs/rnaseq.yaml`)  

### Setting Up the Environment  

1. Clone this repository:  

   ```sh
   git clone https://github.com/anamikagiri/RNAseq.git
   cd rnaseq

2. Activate the snakemake environment:
  
   mamba activate snakemake

3. Usage

   Running the Pipeline:
   For cluster execution (e.g., SLURM):

   snakemake --profile slurm_prof

4. Output

   The pipeline generates:

   Trimmed Reads: results/trimmed/ 

   Aligned Reads: results/alignment/ 
 
   Trancripts Counts : results/salmon/

   Gene Counts: results/featurecounts/ 

   QC Reports: qc/
