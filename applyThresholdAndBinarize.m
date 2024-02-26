function W_thr_b = applyThresholdAndBinarize(Z, thr)
% applyThresholdAndBinarize Applies a proportional threshold and binarizes the connectivity matrices.
%
% Syntax: W_thr_b = applyThresholdAndBinarize(Z, thr)
%
% Inputs:
%   Z - A 3D matrix representing the connectivity data, with dimensions [numNodes x numNodes x nParticipants],
%       where 'numNodes' is the number of nodes in the graph and 'nParticipants' is the number of participants.
%   thr - A vector of threshold values to be applied to the connectivity matrices.
%
% Outputs:
%   W_thr_b - A 4D matrix containing the thresholded and binarized connectivity data for each participant
%             and each threshold value, with dimensions [numNodes x numNodes x nParticipants x length(thr)].

% Number of participants
n = size(Z, 3);

% Initialize the output matrix
W_thr_b = zeros(size(Z, 1), size(Z, 2), n, length(thr));

% Loop over threshold values and participants
for itThr = 1:length(thr)
    for itPp = 1:n
        % Threshold and binarize
        W_thr = threshold_proportional(Z(:,:,itPp), thr(itThr));
        W_thr_b(:,:,itPp,itThr) = W_thr > 0; % Binarize
    end
end

% Convert matrix from logical to double
W_thr_b = double(W_thr_b);
end
