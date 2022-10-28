# Alcoholic Hepatatits Liver Genome-Scale Metabolic Model Study

This repository contains 9 files.
- metabolic-Tasks_Essential.xlsx: Excel file containing the metabolic task list needed for generating the GEMs
- GenerateDiseaseGEM.m: MATLAB script that generates all GEMs using the tINIT algorithm
- StructuralComparisons.m: MATLAB script that calculates the structural comparisons between all GEMs using the compareMultipleModels function
- RunSpotEflux2.m: MATLAB script that calculates the model objective functions and performs flux balance analysis
- GEM_MATLABfunctions.zip: Zipped file containing all MATLAB functions and scripts required for running the code in files GenerateDiseaseGEM.m, StructuralComparisons.m, RunSpotEflux2.m
- GEMs.zip: Zipped file containing all 8 GEM files in .mat format
- FigureReplication.R: R script for replicating the Figures in the present manuscript
- FigureReplicationFiles.zip: Zipped file containing all necessary files for Figure replication
- TenSimpleRules.docx: File containing the document describing the self-assessed conformance of our models to the state of the art Ten Simple Rules for Credible Practice in Modeling and Simulation in Healthcare

# Dependencies
- MATLAB (ver. R2021a and higher)
- RStudio (ver. 4.2.1 and higher)
- [Gurobi Optimizer](https://www.gurobi.com/downloads/gurobi-optimizer-eula/) for MATLAB (ver. 8.1.1)

# MATLAB Toolboxes

# Rtudio packages
- dataVisEasy – 0.3.2
- limma – 3.52.4
- ica - 1.0-3
- matrixStats - 0.62.0
- pheatmap - 1.0.12
- pcaMethods - 1.88.0
- tidyverse - 1.3.2
- dplyr - 1.0.10
- reshape2 - 1.4.4
- stats - 4.1.1
- ape – 5.6-2
