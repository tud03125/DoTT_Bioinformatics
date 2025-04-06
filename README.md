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
Path to the GTF annotation file (e.g., for mm39: ```/path/to/mm39_RefSeq.gtf```; for hg38: ```/path/to/hg38.knownGene.gtf```).

**--bam-files**
A list of BAM file paths (space-separated). The order of the files must correspond to the sample conditions provided with ```--conditions```.

**--species**
Species option. For mouse use ```mm39```; for human use ```hg38``` or ```hg19```.

**--extension**
Fixed extension length in bases (e.g., ```10000```).

**--output-dir**
Directory where all output files will be written.

**--conditions**
Comma-separated list of condition labels for each BAM file.

   For simulated (mouse) data, for example: ```Fasted,Fasted,Fasted,Fasted,Fasted,HCD,HCD,HCD,HCD,HCD```

   For human data, for example: ```mock,mock,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1```

## Optional Arguments
**--dynamic**
Enable dynamic region extension (if not provided, a fixed extension is used).

**--kgx-file**
Path to the kgXref mapping file (used for human GTFs).

**--run_gsea**
Enable generation of a GSEA pre-ranked list from the DESeq2 results.

**--bootstrap, --n_boot, --consensus_threshold**
Enable bootstrapping in the DESeq2 analysis.

   ```--bootstrap``` is a flag (include it to enable bootstrapping).
   
   ```--n_boot``` sets the number of bootstrap iterations (default is 100).
   
   ```--consensus_threshold``` is the fraction (default is 0.5) required for a gene to be considered consensus.

**--unsupervised_ml**
Enable the unsupervised ML analysis module, which uses the DESeq2 results along with the provided conditions to assess replicate consistency and perform enrichment comparisons.

**--supervised_ml**
Enable supervised ML analysis (which compares DESeq2 results to simulation ground truth and trains an ML classifier).

**--sim_tx_info**
Path to the simulation ground truth file (e.g., ```sim_tx_info.txt```).
**Required if ```--supervised_ml``` is used.**

**--experimental_condition**
The label for the experimental condition. This value is used by ML modules to separate experimental vs. control samples.
**Required if** ```--supervised_ml``` **or** ```--unsupervised_ml``` **is used.**

## Example 1: Simulated Mouse (mm39) Test with Bootstrapping, GSEA, and Supervised ML

```
cd /path/to/DoTT_Bioinformatics
python3 main.py \
  --gtf-file /path/to/mm39_RefSeq.gtf \
  --bam-files /path/to/simulated_reads/STAR_sample_01_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_02_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_03_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_04_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_05_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_06_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_07_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_08_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_09_Aligned.sortedByCoord.out.bam \
              /path/to/simulated_reads/STAR_sample_10_Aligned.sortedByCoord.out.bam \
  --species mm39 \
  --extension 10000 \
  --output-dir DoTT_simulation_test \
  --conditions Fasted,Fasted,Fasted,Fasted,Fasted,HCD,HCD,HCD,HCD,HCD \
  --bootstrap \
  --n_boot 100 \
  --consensus_threshold 0.5 \
  --run_gsea \
  --supervised_ml \
  --experimental_condition HCD \
  --sim_tx_info /path/to/simulated_reads/sim_tx_info.txt
```

## Example 2: Human Data (hg38) Test with Bootstrapping and GSEA

```
cd /path/to/DoTT_Bioinformatics
python3 main.py \
  --gtf-file /path/to/hg38.knownGene.gtf \
  --bam-files /path/to/human/dataset/Total_RNA_mock/SRR1523653_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_mock/SRR1523667_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523654_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523655_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523656_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523657_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523668_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523669_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523670_Aligned.sortedByCoord.out.bam \
              /path/to/human/dataset/Total_RNA_Herpes_simplex_virus_1_strain_17/SRR1523671_Aligned.sortedByCoord.out.bam \
  --species hg38 \
  --extension 10000 \
  --output-dir DoTT_HSV-1_mock_test \
  --conditions mock,mock,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1,HSV-1 \
  --bootstrap \
  --n_boot 100 \
  --consensus_threshold 0.5 \
  --run_gsea
```
