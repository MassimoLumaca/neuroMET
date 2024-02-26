function [Local_Efficiency_roi, Clustering_Coefficient_roi, Betweenness_Centrality_roi, E_global_roi] = computeGraphTheoryMetrics(thr, Z, W_thr_b)
% computeGraphTheoryMetrics Computes graph theory metrics for given thresholds and connectivity matrices.
%
% Syntax:
% [Local_Efficiency_roi, Clustering_Coefficient_roi, Betweenness_Centrality_roi, E_global_roi] = computeGraphTheoryMetrics(thr, Z, W, W_thr_b)
%
% Inputs:
%   thr - A vector of threshold values.
%   Z - Original connectivity matrices (unused in the function but included for completeness).
%   W_thr_b - Binarized thresholded connectivity matrices.
%
% Outputs:
%   Local_Efficiency_roi - Local efficiency for each ROI.
%   Clustering_Coefficient_roi - Clustering coefficient for each ROI.
%   Betweenness_Centrality_roi - Betweenness centrality for each ROI.
%   E_global_roi - Global efficiency for each ROI.

% Initialize variables
numpp = size(Z, 3); % Number of participants
totalIterations = size(thr,2) * numpp;
N = size(Z,1); % Number of n

% Compute graph theory metrics
h = waitbar(0, 'Processing...'); % Initialize waitbar
for itThr = 1:length(thr)
    disp(['Compute graph theory metrics through threshold ' num2str(thr(itThr))]);
    for it = 1:numpp
        A = W_thr_b(:,:,it,itThr); % binary matrix

        % Network-level measures with BCT functions
        
        % Global efficiency
        % Calculate shortest path length of the binary matrix
        D = distance_bin(A);
        % Calculate the inverse of D, avoiding division by zero
        iD = zeros(size(D));
        iD(D > 0) = 1 ./ D(D > 0);
        % Calculate the global efficiency for each node
        for n = 1:N
            E_global(n) = sum(iD(n, :)) / (N - 1);
        end
        [E_global_roi(it,:,itThr)] = E_global;
        
        Local_Efficiency_roi(it,:,itThr) = efficiency_bin(A, 1); % Local efficiency
        Clustering_Coefficient_roi(it,:,itThr) = clustering_coef_bu(A); % Clustering coefficient
        Betweenness_Centrality_roi(it,:,itThr) = betweenness_bin(A'); % Betweenness Centrality

        % Update progress on waitbar
        progress = ((itThr-1) * numpp + it) / totalIterations;
        waitbar(progress, h, sprintf('Processing %d%%', floor(progress * 100)));
    end
end

% Close the waitbar
close(h);
end

