function Z = selectNodesForGraph(Z, nodes)
% selectNodesForGraph Selects specific nodes from a connectivity matrix for an undirected graph.
%
% Syntax: Z = selectNodesForGraph(Z, nodes)
%
% Inputs:
%   Z - A 3D matrix representing the connectivity data, with dimensions [numnodes x numnodes x n],
%       where 'numnodes' is the number of nodes in the original graph and 'n' is the number of participants or time points.
%   nodes - A vector containing the indices of nodes to be selected from the graph.
%
% Outputs:
%   Z - A 3D matrix (Znerd) containing the connectivity data for the selected nodes,
%       with dimensions [length(nodes) x length(nodes) x n].

% Number of nodes in the selection and the number of matrices (participants or time points)
numSelectedNodes = length(nodes);
n = size(Z, 3);

% Initialize the output matrix with zeros
Znerd = zeros(numSelectedNodes, numSelectedNodes, n);

% Populate the new matrix with connectivity data for the selected nodes
for it = 1:n % Iterating through each participant or time point
    n1_loop = 0;
    for n1 = nodes % Iterating through the selected nodes
        n1_loop = n1_loop + 1;
        n2_loop = 0;
        for n2 = nodes
            n2_loop = n2_loop + 1;
            Znerd(n1_loop, n2_loop, it) = Z(n1, n2, it);
        end
    end
end

% Return the new matrix with selected nodes' connectivity data
Z = Znerd;
end

