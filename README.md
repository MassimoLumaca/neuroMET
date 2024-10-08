# Frontoparietal network topology as a neuromarker of music perceptual abilities

## Overview

This repository contains scripts designed for the complex-network analysis of structural and functional brain-connectivity data and for assessing their relationship to musical abilities (MET) and working memory (WMI). The scripts are extensively commented to ensure the replicability of the graph theory results detailed in our manuscript and its supplementary materials. 
Data cannot be shared publicly as it is part of an ongoing study, and thus considered unanonymized under Danish law. Researchers who wish to access the data may contact Dr Kristian Sandberg (kristian.sandberg@cfin.au.dk) at The Center of Functionally Integrative Neuroscience and/or The Technology Transfer Office (TTO@au.dk) at Aarhus University, Denmark, to establish a data sharing agreement. After permission has been given by the relevant ethics committee, data will be made available to the researchers for replication purposes.

## System requirements

To use these scripts, users must install specific softwares:

- **MATLAB**: A platform for numerical computation, visualization, and programming. Installation instructions for Mac, Windows, and Linux are available at [MathWorks Installation Guide](https://it.mathworks.com/help/install/ug/install-products-with-internet-connection.html)
- **SPM12**: software package designed for the analysis of brain imaging data sequences
- **Conn Toolbox**: An open-source, MATLAB-based software for analyzing functional connectivity Magnetic Resonance Imaging (fcMRI). Installation instructions for Mac, Windows, and Linux can be found at [Conn Installation Guide](https://web.conn-toolbox.org/resources/conn-installation)
- **Brain Connectivity Toolbox**: a MATLAB toolbox for complex brain-network analysis. To use, download BCT.zip from [BCTnet](https://sites.google.com/site/bctnet), extract its contents, navigate to these contents in MATLAB (or add the directory to the MATLAB path), and execute functions from the MATLAB command window.
- **R**: A programming language and free software environment for statistical computing and graphics. Download the latest R installer for your operating system (Windows, Mac, or Linux) from the official CRAN mirrors at (https://cran.r-project.org/mirrors.html). Run the installer and follow the on-screen instructions to complete the setup.
- **RStudio**: An integrated development environment (IDE) for R. After installing R, download the free RStudio Desktop version from https://posit.co/download/rstudio-desktop/. Run the installer for your platform and follow the installation wizard, keeping the default options. 

Our analysis was tested on Matlab 2016b, Conn toolbox version 2021a, the Brain Connectivity Toolbox (version 2019-03-03), and RStudio (version 2024.04.0+735).

## Installation guide and demo

1. Download this repository as a ZIP file and unzip it into the desired directory.
2. Add the directory path to MATLAB.
3. Set the MATLAB current directory to the repository folder. It contains four main scripts:
	- ``functional_MET.m``: Performs graph theory analysis of functional connectivity and multiple regression statistics with MET or Emotion, controlling for age, sex, and musical training.
	- ``structural_MET.m``: Performs graph theory analysis of structural connectivity and multiple regression statistics with MET or Emotion, controlling for age, sex, and musical training.
	- ``functional_WMI.m``: Performs graph theory analysis of functional connectivity and multiple regression statistics with WMI, controlling for age, sex, and musical training.
	- ``structural_WMI.m``: Performs graph theory analysis of structural connectivity and multiple regression statistics with WMI, controlling for age, sex, and musical training.
	
In the code snippets below, replace ‘name_of_the_script.m’ with the desired script name (e.g., ‘structural_WMI.m’) and replace ‘network_of_interest’ with either ‘frontoparietal’ (for analysis on the frontoparietal network) or ‘occipital’ (for analysis on the occipital network). In ``functional_MET.m`` and ``structural_MET.m`` is possible to change the behavioural variable of interest (e.g., 'PercentageTotalMETScore' from the MET, or 'F5' from the Gold-MSI).

4. Run the following commands in the MATLAB command line:
```
currentDir = pwd;
graph_main = fullfile(currentDir,’name_of_the_script.m’);
target_net = {‘network_of_interest’};
run(graph_main)
```

The execution will generate a table of F-values in the command window (Suppl. Tables 5-12) and two variables in the workspace:

- ``results``: Includes computed adjacency matrix, threshold range, selected ROIs, individual behavioural measures, and graph metrics.
- ``stats``: Contains final statistical results. Refer to stats.p_two_side_fdr for main manuscript results.

### Interpretation
To map p-values in ``stats.p_two_side_fdr`` to the correct brain regions, match the p-value index with the same index in ``results.rois``. For example, if the p-value of interest is located in the second column of ``stats.p_two_side_fdr``, the associated brain region will be listed in the second row of ``results.rois``. 

## Other scripts

The repository also contains auxiliary scripts for multivariate outlier analysis with PCA (``multivariate_outlier_detection`` folder) and mediation analysis (``mediation_analysis`` folder), each thoroughly commented for ease of use. 
