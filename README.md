# Intro to scRNAseq:
Introduction to single cell RNAseq data analysis and interpretation coursework material provided by Harvard Informatics Team.

#### Prepare course material
Get the material in this repository
```
git clone https://github.com/yasinkaymaz/singleCellRNAseq_workshop.git
```
Prepare R/Rstudio environment

For Seurat, follow the instructions: https://satijalab.org/seurat/install.html
or simply

```
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')

source("https://bioconductor.org/biocLite.R")
biocLite("SIMLR")

library(Seurat)
library(Matrix)
library(dplyr)
library(SIMLR)
```

Download the datasets into your ~/Downloads folder.

```
cd ~/Downloads
wget https://ifx.rc.fas.harvard.edu/pub/bionano2018/scRNAseqWS.zip
unzip scRNAseqWS.zip
```

Follow the tutorial in workshop-workflow.html
