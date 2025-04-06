# DoTT_Bioinformatics
A bioinformatics pipeline for analyzing disruption of transcriptional termination (DoTT) using Python and R (via rpy2), designed for both mouse and human datasets.

## Overview
This repository contains Python pipeline for performing DoTT (disruption of transcriptional termination) analysis. It integrates multiple processing steps such as generating SAF files from GTF annotations, running featureCounts, filtering counts, differential expression analysis with DESeq2 (via rpy2), and merging significant gene regions. The pipeline supports both mouse (mm39) and human (hg38/hg19) datasets with flexible input parameters.

## Features
- **Species-specific processing:** Options for mm39 (mouse) and hg38/hg19 (human) with different gene ID strategies.
- **Flexible region extension:** Choose between fixed and dynamic extension based on read density.
- **Differential expression analysis:** Runs DESeq2 via rpy2 and creates MA and Volcano plots.
- **Modular design:** Easily configurable via command-line arguments.

## Installation
1. **Clone the Repository:**
   ```bash
   gh repo clone tud03125/DoTT_Bioinformatics
   cd DoTT_Bioinformatics
   
2. **Pre-requisites (using pip):**
   ```
   pip install -r requirements.txt
   ```
   Or

   **Pre-requisites (using conda):**
   ```
   conda env create -f environment.yml
   ```

   **R Package Installation**

   This pipeline requires the following R packages:
      -DESeq2 (install via Bioconductor)
      -EnhancedVolcano (install via Bioconductor)
      -argparse, ggplot2, ggrepel (install via CRAN)
   
   To install these packages, you can run the following in R:

   ```
   # Install CRAN packages:
   install.packages(c("argparse", "ggplot2", "ggrepel"), repos = "https://cran.rstudio.com")
   
   # Install Bioconductor packages:
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("DESeq2", "EnhancedVolcano"))
   ```

4. **Usage**

## Required Arguments
**--gtf-file**
Path to the GTF annotation file (e.g., for mm39: /home/tud03125/pipeline/mm39_RefSeq.gtf; for hg38: /home/tud03125/pipeline/hg38.knownGene.gtf).

**--bam-files**
A list of BAM file paths (space-separated). The order of the files must correspond to the sample conditions provided with --conditions.

**--species**
Species option. For mouse use mm39; for human use hg38 or hg19.

**--extension**
Fixed extension length in bases (e.g., 10000).

**--output-dir**
Directory where all output files will be written.

**--conditions**
Comma-separated list of condition labels for each BAM file.

For simulated (mouse) data, for example: Fasted,Fasted,Fasted,Fasted,Fasted,HCD,HCD,HCD,HCD,HCD

For human data, for example: mock,mock,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1

## Optional Arguments
**--dynamic**
Enable dynamic region extension (if not provided, a fixed extension is used).

**--kgx-file**
Path to the kgXref mapping file (used for human GTFs).

**--run_gsea**
Enable generation of a GSEA pre-ranked list from the DESeq2 results.

**--bootstrap, --n_boot, --consensus_threshold**
Enable bootstrapping in the DESeq2 analysis.

--bootstrap is a flag (include it to enable bootstrapping).

--n_boot sets the number of bootstrap iterations (default is 100).

--consensus_threshold is the fraction (default is 0.5) required for a gene to be considered consensus.

**--supervised_ml**
Enable supervised ML analysis (which compares DESeq2 results to simulation ground truth and trains an ML classifier).

**--sim_tx_info**
Path to the simulation ground truth file (e.g., sim_tx_info.txt).
Required if --supervised_ml is used.

**--experimental_condition**
The label for the experimental condition. This value is used by ML modules to separate experimental vs. control samples.
Required if --supervised_ml is used.
