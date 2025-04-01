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

3. **Usage**

To run the DoTT Bioinformatics Pipeline, use the following command-line template:

 ```
python DOTT_Bioinformatics_full_pipeline_GitHub.py \
  --gtf-file [path/to/gtf_file] \
  --bam-files [path/to/bam_file1] [path/to/bam_file2] ... [path/to/last_bam_file] \
  --species [mm39 | hg38 | hg19] \
  --extension [extension_length_in_bases] \
  --output-dir [path/to/output_directory] \
  --conditions "[condition1],[condition2],..., [conditionN]"
 ```
