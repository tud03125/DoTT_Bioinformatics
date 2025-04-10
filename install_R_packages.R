# For CRAN packages
required_cran_packages <- c("argparse", "ggplot2", "ggrepel")
new_cran_packages <- required_cran_packages[!(required_cran_packages %in% installed.packages()[, "Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages, repos="https://cran.rstudio.com")

# For Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
required_bioc_packages <- c("DESeq2", "EnhancedVolcano")
new_bioc_packages <- required_bioc_packages[!(required_bioc_packages %in% installed.packages()[, "Package"])]
if(length(new_bioc_packages)) BiocManager::install(new_bioc_packages)
