# ccRCC Molecular Subtyping Study

## Setup

Before running the code, make sure to install all required R packages. Specifically, the "MOVICS" package is vital. Install "MOVICS" from its GitHub repository:

MOVICS installation
devtools::install_github("xlucpu/MOVICS")

graphql
Copy code

## Data and Folder Structure

After setting up the packages, define the paths for tumor data, input data, figures, and scripts:

```R
tumor.path <- "set your own path"
setwd(tumor.path) #create directory
The code will automatically create the following directories for you:

InputData
Figure1, Figure2, Figure3, Figure5, Figure6, FigureS
Tables
Scripts
Please ensure the tumor.path directory exists or the code will create it for you.

Dependencies
Load the necessary R libraries:

R
Copy code
library(MOVICS)
library(aplot)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(wesanderson)
Analysis and Figures
Figure S1: This part deals with clustering number of TCGA-KIRC and consensus heatmap. Sample code:
R
Copy code
load("./InputData/TCGA_KIRC.fpkm.rda")
load("./InputData/mo.data.rda")
#... [rest of the code for Figure S1]
Figure 1: This part deals with comprehensive heatmap of consensusMOIC, PFI Kaplan-Meier curve, and OS Kaplan-Meier curve.
R
Copy code
#... [code for Figure 1]
Figure S2: Focuses on DEGs and different pathways.
R
Copy code
#... [code for Figure S2]
Figure 2: This part deals with immune gene sets of interest heatmap and immune checkpoints across 3 groups.
R
Copy code
#... [code for Figure 2]
Further Analysis
Continue your analysis as structured in your script. Remember to check each figure and table path to ensure the outputs are saved correctly.

Contribution
If you'd like to contribute or have any issues with the current code, please raise an issue or submit a pull request.

License
[Specify any license or copyright notices]

Contact
[Your contact information]
