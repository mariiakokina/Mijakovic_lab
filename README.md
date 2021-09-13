# Instructions for Working with Day 4 - MSstats Exercises

Before you start working on the exercises, you need to make sure that:

1. R is installed ([https://cloud.r-project.org](https://cloud.r-project.org))
2. RStudio is installed ([https://rstudio.com/products/rstudio/download/](https://rstudio.com/products/rstudio/download/))
3. The following R packages are installed from Bioconductor -- execute the code below in the R console
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSstats")
```


## 

This directory contains scripts and example files from YouTube workshop for MS proteomics data analysis with R package MSstats explaining how to:

1. Read output from Skyline and MaxQuant
2. Perform Data exploration
3. Convert to MSstats required format
4. Data processing (reformat and normalize data)
4. Visualizing data (QC plots)
5. Generate contrast matrix for three comparisons 
6. Inference (subset of significant comparisons; Visualization of differentialy abundant proteins)
7. Make volcano plot
8. Calculating statistical power
9. Visualizing the relationship between desired fold-change and mininum sample size number

Here is the Youtube tutorial:

https://www.youtube.com/watch?v=7chKhQ2UFx4&list=PL2u38g_AG4MH4bc9AZjckRvJhltleGziX

Download the original workshop repository (or clone if you know Git): [https://github.com/ZenBrayn/asms_2020_fall_workshop/archive/main.zip](https://github.com/ZenBrayn/asms_2020_fall_workshop/archive/main.zip)
