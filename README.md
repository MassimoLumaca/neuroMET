# Frontoparietal network topology as a neuromarker of music perceptual abilities

## Overview

This repository contains MATLAB scripts designed for the complex-network analysis of structural and functional brain-connectivity data and for assessing their relationship to musical abilities (MET) and working memory (WMI). The scripts are extensively commented to ensure the replicability of the graph theory results detailed in our manuscript and its supplementary materials. 

## System requirements

To use these scripts, users must install specific softwares:

- **MATLAB**: A platform for numerical computation, visualization, and programming. Installation instructions for Mac, Windows, and Linux are available at [MathWorks Installation Guide](https://it.mathworks.com/help/install/ug/install-products-with-internet-connection.html)
- **Conn Toolbox**: An open-source, MATLAB-based software for analyzing functional connectivity Magnetic Resonance Imaging (fcMRI). Installation instructions for Mac, Windows, and Linux can be found at [Conn Installation Guide](https://web.conn-toolbox.org/resources/conn-installation)
- **Brain Connectivity Toolbox**: a MATLAB toolbox for complex brain-network analysis. To use, download BCT.zip from [BCTnet](https://sites.google.com/site/bctnet), extract its contents, navigate to these contents in MATLAB (or add the directory to the MATLAB path), and execute functions from the MATLAB command window.

Our analysis was tested on Matlab 2016b, Conn toolbox version 2021a, and the Brain Connectivity Toolbox (version 2019-03-03).

## Installation guide and demo

1. Download this repository as a ZIP file and unzip it into the desired directory.
2. Add the directory path to MATLAB.
3. Set the MATLAB current directory to the repository folder. It contains four main scripts:
	- ``functional_MET.m``: Performs graph theory analysis of functional connectivity and multiple regression statistics with MET, controlling for age, sex, and musical training.
	- ``structural_MET.m``: Performs graph theory analysis of structural connectivity and multiple regression statistics with MET, controlling for age, sex, and musical training.
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
