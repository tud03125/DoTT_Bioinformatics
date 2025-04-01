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
   
2. **Pre-requisites**
   ```bash
   pip install -r requirements.txt

   conda env create -f environment.yml
