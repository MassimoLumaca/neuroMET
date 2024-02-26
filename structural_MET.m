% -------------------------------------------------------------------------------------------------------------------
% MATLAB SCRIPT FOR GRAPH THEORY ANALYSIS IN STRUCTURAL CONNECTOMES AND RELATIONSHIP TO MUSICAL EAR TEST
% -------------------------------------------------------------------------------------------------------------------
%
% Developed by: Massimo Lumaca (2023), Associate Professor, Center For
% Music in the Brain, Department of Clinical Medicine, Aarhus University
%
% Overview:
% This script is designed to perform graph theory analyses on static connectomes. Its primary aim is to
% investigate the correlations between neural network metrics and the
% percentage of Total MET (Musical Ear Test) scores. 
%
% Key Features:
% - Utilization of graph theory metrics to analyze connectomes.
% - Correlation of neural network metrics with behavioural scores.
% - Inclusion of node-level measures such as Local Efficiency, Clustering Coefficient, Betweenness Centrality, and 
%   Global Efficiency to provide a comprehensive analysis of brain network characteristics.
%
% Dependencies:
% - Brain Connectivity Toolbox (BCT): A comprehensive toolbox for complex network analysis, utilized for calculating 
%   various graph theory metrics.
% - Conn: A functional connectivity analysis software, employed for preprocessing, analysis, and visualization of
%   functional connectome data. Note that results may vary slightly based on the CONN version used.
%
% Usage Notes:
% - Ensure that both BCT and Conn software packages are properly installed and accessible in your MATLAB environment 
%   prior to running this script.
% - Adaptations may be required to align with specific research objectives or to accommodate updates in the
%   dependencies and paths

%% Select the correct parameters
clc % Clear command window

thr =[0.15:0.01:0.30]; % range of proportional threshold applied to connectome matrices
behav={'PercentageTotalMETScore'}; % Define behavior variable for correlation analysis

% Participant removal option
removeparticipant={'yes'};
% Set 'yes' if participants need to be removed according to multivariate outlier detection.
% Otherwise, set to 'no'.


%% Selection of the network of interest

% Determine which nodes to include based on the selected network
if strcmp(target_net,'frontoparietal')
    % Frontoparietal network nodes
    nodes=[15 16 25 26 27 53 54 56 89 90 99 100 101 127 128 130];
elseif strcmp(target_net,'occipital')
    % Occipital network nodes
    nodes= [2 19 20 21 22 42 59 65 76 93 94 95 96 116 133 139]; 
end
%% Remove participants

% Conditionally remove participants if specified
if strcmp(removeparticipant,'yes')
    % IDs of participants to remove based on outlier detection
    rmPp = {'0052', '0078','0104', '0130', '0135', '0143', '0149', '0164', '0165', '0196', '0229', '0238', '0260', '0266', '0278', '0283'};   
end

% Participants whose structural connectome construction failed: 0052 0104
% 0165 0229 0260 0266 0278. The others are those that result as outliers in
% the multivariate outlier detection analysis.

%% Load connectomes
% Load functional connectivity matrices and relevant node information

% Determine the directory of the current script
currentScriptDir = fileparts(mfilename('fullpath'));

% Construct paths to necessary directories
connectomesDir = fullfile(currentScriptDir, 'connectomes');
% Paths to specific files within the directories
connectomeFilePath = fullfile(connectomesDir, 'structural_connectome_233.mat');
% Load the specified files
load(connectomeFilePath)

%% Remove connectome matrix of participants not included in the analysis

% Get the list of participant names to be removed
names = rmPp; % 'rmPp' should contain the names of participants to be removed

% Load the file containing folders for participant graphs
load('folders_pp_graphs.mat') % This loads a variable of the participants' folders

% Initialise an empty array to store indices of participants to remove
idx = [];

% Loop through each participant in the directory to identify ones to remove
for k = 1:numel(images_dir) % Iterate over the array of participant directories
    if any(strcmp(names, images_dir(k).name)) % Check if the current directory's name is in the removal list
        idx = [idx; k]; % If so, add the index to the list of indices to remove
    end
end

% Remove the specified slices from the connectome matrix
matrix(:, :, idx) = []; % Removes the matrices corresponding to participants in 'idx'
%% GRAPH THEORY
%% extract name and number nodes Destrieux parcellation

% Construct the path to the 'NodeNames_and_coordinates' directory from the current script's location
NodeNamesAndCoordinatesDir = fullfile(currentScriptDir, 'NodeNames_and_coordinates');
% Construct the path to your specific node names and coordinates files
NodeNamesPath = fullfile(NodeNamesAndCoordinatesDir, 'Destrieux.coordinates.mat');
CoordinatesPath = fullfile(NodeNamesAndCoordinatesDir, 'DestrieuxCorticalNodes148.mat');
load(NodeNamesPath);
load(CoordinatesPath);

% sort myData by a specific column so to have ROI names corresponding with
% those of Conn
myData = sortrows(myData,7);

% Names of ROIs with the same order used in CONN
NodeNames = myData(:,1);
node148Destrieux = NodeNames;


% extract relevant ROI names (Destrieux)
NodeNames = NodeNames(nodes);
numnodes=length(nodes);

%% Node selection for an undirected graph
disp(['Node selection']);
% Select nodes for graph analysis based on the defined network nodes
Z = selectNodesForGraph(matrix, nodes);
%% Thresholding and binarization
disp(['Thresholding and binarization']);
% Apply threshold and binarize the matrix Z using specified thresholds
W_thr_b = applyThresholdAndBinarize(Z, thr);
%% Compute graph theory metrics

[Local_Efficiency_roi, Clustering_Coefficient_roi, Betweenness_Centrality_roi, E_global_roi] = computeGraphTheoryMetrics(thr, Z, W_thr_b);

%% Use NaN for nodes with less than 2 edges

% Initialize number of participants based on the third dimension of Z
n = size(Z, 3);
% Iterate over thresholds, participants, and nodes to sum connections
for itThr=1:size(thr,2)
for pp=1:n
    for nodesK=1:numnodes 
        t(pp,nodesK,itThr)=sum(W_thr_b(nodesK,:,pp,itThr),2); % Sum of edges for each node/participant at each threshold       
    end
end
end

% Assign NaN to nodes with less than 2 edges
NaNmath_1 = t<2; 
Local_Efficiency_roi(NaNmath_1) = NaN;
Clustering_Coefficient_roi(NaNmath_1) = NaN;
%% Average across thresholds
 disp(['Average across thresholds'])
 % If there are multiple thresholds, calculate the mean of graph metrics across thresholds
if size(thr,2) > 1
     Local_Efficiency_roi = mean(Local_Efficiency_roi,3);
     Clustering_Coefficient_roi = mean(Clustering_Coefficient_roi,3);
     Betweenness_Centrality_roi = mean(Betweenness_Centrality_roi,3);
     E_global_roi = mean(E_global_roi,3);
end
%% % Store results in a structured format
results.Z = Z;
results.Z_thr = W_thr_b;
results.thr=thr;
results.rois=NodeNames;
results.measures.Local_Efficiency_roi=Local_Efficiency_roi;
results.measures.Clustering_Coefficient_roi=Clustering_Coefficient_roi;
results.measures.Betweenness_Centrality_roi=Betweenness_Centrality_roi;
results.measures.E_global_roi=E_global_roi;
%% Load behavioural variables
disp(['Load variables']);

% Define the original number of participants
N=233;

behaviourDir = fullfile(currentScriptDir, 'behavioural_data');
behaviourFilePath = fullfile(behaviourDir, 'behavior_update.xlsx');
behavior_table=readtable(behaviourFilePath);

% load mappings brain-behaviour
load('Index_Var_structural.mat')

%% Setup Names Vars in Batch
results.covariates.effect_names ={'AllSubjects','Age', 'Gender', 'MelodyScore','PercentageMelodyScore','RhythmScore','PercentageOfRhythmScore',...
    'TotalMETScore','PercentageOfTotalMETScore','F1','F2','F3','F4','F5','GS'};

% Load behavioral data from table
results.covariates.effect{1}=ones(N,1); %allSubjects

% Extract and assign other covariates based on effect names
for ii = 2:size(results.covariates.effect_names,2)
 namevar = results.covariates.effect_names{1,ii};
 column_values = behavior_table.(namevar)(IndexVar);
 results.covariates.effect{ii}=column_values;
end
%% Remove Participants if Needed


% Remove participants
if strcmp(removeparticipant,'yes')
    for i = 1:size(results.covariates.effect,2)
        CurrCell = results.covariates.effect{1,i};
        
        % remove the specific rows
        CurrCell(idx,:) = [];
        
        % update the cell structure
        results.covariates.effect{1,i} = CurrCell;
    end
end

%% Statistical Analysis Preparation
disp(['Statistical analysis']);

%take fieldnames
fns=fieldnames(results.measures);

% Prepare for GLM analysis by setting up covariates and contrasts
if strcmp(behav,'PercentageTotalMETScore')
    
    % Set the GLM with specific covariates (e.g., intercept, age, gender, PercentageTotalMETScore, Musical Training)
    X = [results.covariates.effect{1,1} results.covariates.effect{1,2} results.covariates.effect{1,3} results.covariates.effect{1,9} results.covariates.effect{1,12}];
    
    C = [0 0 0 1 0]; % contrast vector: intercept, age, sex, GraphMetrics, MusicalTraining
    M = 1;
    D = [0 0 0 0 0];
    
    disp(['Effect of PercentageTotalMETScore, while correcting for age, gender, and training']);
end
%% Remove NaN Participants for Nodal-Level Measures

for ll=1:numnodes
    
% Process each nodal-level measure to remove NaN values, ensuring clean data for statistical analysis
% Local Efficiency
meas_temp_locEff=results.measures.Local_Efficiency_roi(:,ll);
meas_temp_locEff_Nan= meas_temp_locEff(~isnan(meas_temp_locEff));
results.measures.measuresNaN.Local_Efficiency_roi.scores{:,ll}=meas_temp_locEff_Nan;
 
meas_temp_ClustCoeff=results.measures.Clustering_Coefficient_roi(:,ll);
meas_temp_ClustCoeff_Nan= meas_temp_ClustCoeff(~isnan(meas_temp_ClustCoeff));
results.measures.measuresNaN.Clustering_Coefficient_roi.scores{:,ll}=meas_temp_ClustCoeff_Nan;

meas_temp_BetweennCentral=results.measures.Betweenness_Centrality_roi(:,ll);
meas_temp_BetweennCentral_Nan= meas_temp_BetweennCentral(~isnan(meas_temp_BetweennCentral));
results.measures.measuresNaN.Betweenness_Centrality_roi.scores{:,ll}=meas_temp_BetweennCentral_Nan;

meas_temp_E_global=results.measures.E_global_roi(:,ll);
meas_temp_E_global_Nan= meas_temp_E_global(~isnan(meas_temp_E_global));
results.measures.measuresNaN.E_global_roi.scores{:,ll}=meas_temp_E_global_Nan;

clear meas_temp_locEff meas_temp_locEff_Nan meas_temp_ClustCoeff_Nan meas_temp_BetweennCentral meas_temp_BetweennCentral_Nan
end
%% Run GLM for the graph metrics of interest
% h = effect size (beta); F = test used; p = one-sided p-value; dof = degrees of freedom

% Node-level measures analysis
% Local efficiency
    for itr=1:numnodes % Loop across all nodes
        % Identify NaN values for each participant in the Local Efficiency measure
        NaNpp = isnan(results.measures.Local_Efficiency_roi(:,itr));
        Xmat = X(~NaNpp,:); % Remove participants with NaN values from the design matrix X
        Y = [results.measures.measuresNaN.Local_Efficiency_roi.scores{1,itr}]; % Remove NaN values from the dependent variable 
        [h,F,p,dof] = conn_glm(Xmat,Y,C,M,D); % Perform GLM
        p_twoside_unc = 2*min(p,1-p); % Calculate two-sided uncorrected p-value
        
        % Store the statistics in the stats structure
        stats.(fns{1}).h(itr)=h(:,:,1);
        stats.(fns{1}).F(itr)=F(:,:,1);
        stats.(fns{1}).p(itr)=p(:,:,1);
        stats.(fns{1}).p_twoside_unc(itr)=p_twoside_unc(:,:,1);
        stats.(fns{1}).p_twoside_fdr=conn_fdr(stats.Local_Efficiency_roi.p_twoside_unc)
        stats.(fns{1}).dof(itr)=dof;
        stats.(fns{1}).statsname='T'; % Type of statistical test
        clear Xmat Y h F p dof % Clear temporary variables
    end
    
% Clustering Coefficient
    for itr=1:numnodes 
        
        NaNpp = isnan(results.measures.Clustering_Coefficient_roi(:,itr));
        Xmat = X(~NaNpp,:); 
        Y = [results.measures.measuresNaN.Clustering_Coefficient_roi.scores{1,itr}]; %remove Nan pp from Y
        [h,F,p,dof] = conn_glm(Xmat,Y,C,M,D);
        p_twoside_unc = 2*min(p,1-p); 
        
        stats.(fns{2}).h(itr)=h(:,:,1);
        stats.(fns{2}).F(itr)=F(:,:,1);
        stats.(fns{2}).p(itr)=p(:,:,1);
        stats.(fns{2}).p_twoside_unc(itr)=p_twoside_unc(:,:,1);
        stats.(fns{2}).p_twoside_fdr=conn_fdr(stats.Clustering_Coefficient_roi.p_twoside_unc)
        stats.(fns{2}).dof(itr)=dof;
        stats.(fns{2}).statsname='T';
        clear Xmat Y h F p dof
    end
    
    % Betweenness centrality
    for itr=1:numnodes 
        
        NaNpp = isnan(results.measures.Betweenness_Centrality_roi(:,itr));
        Xmat = X(~NaNpp,:); 
        Y = [results.measures.measuresNaN.Betweenness_Centrality_roi.scores{1,itr}]; %remove Nan pp from Y
        [h,F,p,dof] = conn_glm(Xmat,Y,C,M,D);
        p_twoside_unc = 2*min(p,1-p);
        
        stats.(fns{3}).h(itr)=h(:,:,1);
        stats.(fns{3}).F(itr)=F(:,:,1);
        stats.(fns{3}).p(itr)=p(:,:,1);
        stats.(fns{3}).p_twoside_unc(itr)=p_twoside_unc(:,:,1);
        stats.(fns{3}).p_twoside_fdr=conn_fdr(stats.Betweenness_Centrality_roi.p_twoside_unc)
        stats.(fns{3}).dof(itr)=dof;
        stats.(fns{3}).statsname='T';
        clear Xmat Y h F p dof
    end
   
    
    % E_global ROI
    
    for itr=1:numnodes 
        
        NaNpp = isnan(results.measures.E_global_roi(:,itr));
        Xmat = X(~NaNpp,:); 
        Y = [results.measures.measuresNaN.E_global_roi.scores{1,itr}]; 
        [h,F,p,dof] = conn_glm(Xmat,Y,C,M,D);
        p_twoside_unc = 2*min(p,1-p); 
        
        stats.(fns{4}).h(itr)=h(:,:,1);
        stats.(fns{4}).F(itr)=F(:,:,1);
        stats.(fns{4}).p(itr)=p(:,:,1);
        stats.(fns{4}).p_twoside_unc(itr)=p_twoside_unc(:,:,1);
        stats.(fns{4}).p_twoside_fdr=conn_fdr(stats.E_global_roi.p_twoside_unc)
        stats.(fns{4}).dof(itr)=dof;
        stats.(fns{4}).statsname='T';
        clear Xmat Y h F p dof
    end
clearvars -except results stats NodeNames

%% Create a table with F-values for the results, to be used in creating a heatmap (Supplementary Tables 5-6)

% Extract F-values from the stats structure for each measure
Local_Efficiency = stats.Local_Efficiency_roi.F';
Clustering_Coefficient = stats.Clustering_Coefficient_roi.F';
Betweenness_Centrality = stats.Betweenness_Centrality_roi.F';
E_global = stats.E_global_roi.F';

% Create a table with the extracted F-values, using NodeNames as row names
T = table(Local_Efficiency,Clustering_Coefficient,Betweenness_Centrality,E_global,'RowNames',NodeNames)
disp(T) % Display the table

% Add RowNames as a new variable in the table
T.NodeNames = NodeNames

% Move the 'NodeNames' variable to be the first column in the table
T = [T(:,end), T(:, 1:end-1)];

% Clear all variables except for results and stats for future use
clearvars -except results stats
