# Alcoholic Hepatitis Liver Genome-Scale Metabolic Model Study

This repository contains 9 files.
- **metabolic-Tasks_Essential.xlsx**: Excel file containing the metabolic task list needed for generating the GEMs
- **GenerateDiseaseGEM.m**: MATLAB script that generates all GEMs using the tINIT algorithm
- **StructuralComparisons.m**: MATLAB script that calculates the structural comparisons between all GEMs using the compareMultipleModels function
- **RunSpotEflux2.m**: MATLAB script that calculates the model objective functions and performs flux balance analysis
- **GEM_MATLABfunctions.zip**: Zipped file containing all MATLAB functions and scripts required for running the code in files GenerateDiseaseGEM.m, StructuralComparisons.m, RunSpotEflux2.m
- **GEMs.zip**: Zipped file containing all 8 GEM files in .mat format
- **FigureReplication.R**: R script for replicating the Figures in the present manuscript
- **FigureReplicationFiles.zip**: Zipped file containing all necessary files for Figure replication
- **TenSimpleRules.docx**: File containing the document describing the self-assessed conformance of our models to the state of the art Ten Simple Rules for Credible Practice in Modeling and Simulation in Healthcare

# Dependencies
- MATLAB (ver. R2021a and higher)
- RStudio (ver. 4.2.1 and higher)
- [Gurobi Optimizer](https://www.gurobi.com/downloads/gurobi-optimizer-eula/) for MATLAB (ver. 8.1.1)

# MATLAB Toolboxes
- [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) model (ver. 1.6.0)
- [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox (ver. 2.4.1)

# Rtudio Packages
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

# Usage - Generating GEMS and Performing Flux Balance Analysis
**1. Downloading/Installing toolboxes and organizing workspace**  
**2. Generating genome scale metabolic models (GEMs)**  
**3. Structural Analysis of GEMs**  
**4. Calculating objective functions using SPOT and flux balance analysis using E-Flux2**    
\
Because of the intense computational power necessary for generating the genome scale metabolic models, we suggest running, at minimum, Part 2 - Generating genome scale metabolic models (below) on a server. We utilized a 56-core server with Intel Xeon CPU E5-2697.
\
\
Whether you choose to run this code on a computational server or locally, you must have a working version of MATLAB available. The instructions and example code provided herein are specified for individuals choosing to run the code on a computational server. However, the same principles can be applied for individuals choosing to run the code locally.

**Part 1 - Downloading/Installing toolboxes and organizing workspace**
1. Create a directory named “GEMGeneration”.  
	&emsp;&emsp;&emsp;Example code: mkdir ‘GEMGeneration’
3. Download [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) model (ver. 1.6.0). Place this folder in the “GEMGeneration” directory.  
4. Download [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox (ver. 2.4.1). Place this toolbox in the “GEMGeneration” directory.  
5. Download and install [Gurobi Optimizer](https://www.gurobi.com/downloads/gurobi-optimizer-eula/) (ver. 8.1.1). Link the optimizer to your working version of Matlab using the provided instructions from Gurobi.  
6. Place the following file from FigureReplicationFiles.zip in the “GEMGeneration” directory on the server using a file transfer protocol: “DiseaseData_TPMnorm.txt”.  
	&emsp;&emsp;&emsp;Example code: sftp ‘name of server’  
	&emsp;&emsp;&emsp;put ‘DiseaseData_TPMnorm.txt’ ‘GEMGeneration’  
7. Place the following file in the “GEMGeneration” directory on the server using a file transfer protocol: “GenerateDiseaseGEM.m” and “StructuralComparisons.m”.  
	&emsp;&emsp;&emsp;Example code: sftp ‘name of server’  
	&emsp;&emsp;&emsp;put ‘GenerateDiseaseGEM.m’ ‘GEMGeneration’  
	&emsp;&emsp;&emsp;put ‘StructuralComparisons.m’ ‘GEMGeneration’  
8. Create a directory within the “GEMGeneration” folder named “DiseaseGEM”.  
	&emsp;&emsp;&emsp;Example code: mkdir ‘GEMGeneration/DiseaseGEM’  

**Part 2 - Generating genome scale metabolic models**
1. Since the generation of GEM’s may take a long time, you may want to start a new screen in the server so that the code continues to run even if you disconnect from the server. Learn more about the screen command here: https://linuxize.com/post/how-to-use-linux-screen/   
    &emsp;&emsp;&emsp;Example Code: screen  
2. Start up Matlab on the server  
   &emsp;&emsp;&emsp;Example Code: matlab  
3. Add “GEMGeneration” directory to your Matlab path.  
   &emsp;&emsp;&emsp;Example Code: addpath(genpath(‘GEMGeneration’))  
4. Run ‘GenerateDiseaseGEM.m’ script in Matlab  
    &emsp;&emsp;&emsp;Example Code: GenerateDiseaseGEM  

**Part 3 - Structural Analysis of GEMs**
1. Start up Matlab on the server  
    &emsp;&emsp;&emsp;Example Code: matlab  
2. Add “GEMGeneration” directory to your Matlab path.  
    &emsp;&emsp;&emsp;Example Code: addpath(genpath(‘GEMGeneration’))  
3. Run “StructuralComparisons.m” script in Matlab  
    &emsp;&emsp;&emsp;Example Code: StructuralComparisons  
4. Close Matlab and use a file transfer protocol to put the generated GEM files and structural comparison files in your local computer directory  
   &emsp;&emsp;&emsp;Example Code: quit(matlab)  
   &emsp;&emsp;&emsp;sftp ‘name of server’  
   &emsp;&emsp;&emsp;get ‘GEMGeneration/DiseaseGEM’ .  

**Part 4 – Calculating objective functions using SPOT and flux balance analysis using E-Flux2**
1. Create another directory on the same level as the “GEMGeneration” directory called “FluxBalanceAnalysis”  
     &emsp;&emsp;&emsp;Example code: mkdir ‘FluxBalanceAnalysis’  
2. Place the following file from GEM_MATLABfunctions.zip in the “FluxBalanceAnalysis” directory on the server using a file transfer protocol: “Human-GEM-reactions.txt”  
    Example code: sftp ‘name of server’  
	 &emsp;&emsp;&emsp;put ‘Human-GEM-reactions.txt’ ‘FluxBalanceAnalysis’  
3. Place the following files from GEM_MATLABfunctions.zip in the “FluxBalanceAnalysis” directory on the server using a file transfer protocol: “flux_balance.m”, “solve_milp.m”, “mat_SPOTformat.m”, and “SpotEflux2.m”. and “RunSpotEflux2.m” .  
   	 &emsp;&emsp;&emsp;Example code: sftp ‘name of server’  
	 &emsp;&emsp;&emsp;put ‘flux_balance.m’ ‘FluxBalanceAnalysis’  
 	 &emsp;&emsp;&emsp;put ‘solve_milp.m’ ‘FluxBalanceAnalysis’  
	 &emsp;&emsp;&emsp;put ‘mat_SPOTformat.m’ ‘FluxBalanceAnalysis’  
	 &emsp;&emsp;&emsp;put ‘SpotEflux2.m’ ‘FluxBalanceAnalysis’  
	 &emsp;&emsp;&emsp;put ‘RunSpotEflux2.m’ ‘FluxBalanceAnalysis’  
4. Start up Matlab on the server and add both “GEMGeneration” and “FluxBalanceAnalysis” directories to your Matlab path  
         &emsp;&emsp;&emsp;Example Code: matlab  
	 &emsp;&emsp;&emsp;addpath(genpath(‘GEMGeneration’))  
	 &emsp;&emsp;&emsp;addpath(genpath(‘FluxBalanceAnalysis’))  
5. Run the SPOT and E-Flux2 algorithms by running the ‘RunSpotEflux2.m’ script to calculate the objective function and perform flux balance analysis, respectively. The script saves the generated objective functions, fluxes and models in ‘.mat’ format. It also saves a matrix of all fluxes and a matrix of all fluxes greater than 0.1 to files in ‘.txt’ format.  
     &emsp;&emsp;&emsp;Example Code: RunSpotEflux2  
6. Close Matlab and use a file transfer protocol to put the file containing all fluxes in your local computer directory You may transfer the additional files generated from the RunSpotEflux2 script, however they are not necessary for downstream analysis. All files can be found in the “FluxBalanceAnalysis” directory.  
   &emsp;&emsp;&emsp;Example Code: quit(matlab)  
   &emsp;&emsp;&emsp;sftp ‘name of server’  
   &emsp;&emsp;&emsp;get ‘FluxBalanceAnalysis/DiseaseStateGEM_allfluxes.txt’ .  
Note: The “DiseaseStateGEM_allfluxes.txt” file downloaded from the server is the same as the “DiseaseState_Fluxes.txt” file in FigureReplicationFiles.zip.  

# Usage - Figure Replication
1. Open RStudio
2. Install all necessary packages (make sure the package versions match those from "RStudio Packages")
3. Open FigureReplication.R script and change the working directory to the extracted 'FigureReplicationFiles.zip' folder
4. Run the 'FigureReplication.R' script
