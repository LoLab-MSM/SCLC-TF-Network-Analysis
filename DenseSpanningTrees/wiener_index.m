function [W_ind, D] = wiener_index(T)
    % Reference: K. THILAKAM & A. SUMATHI, HOW TO COMPUTE THE WIENER INDEX OF A
    % GRAPH USING MATLAB, International Journal of Applied Mathematics &
    % Statistical Sciences (IJAMSS), Vol. 2, Issue 5, Nov 2013, 143-148
    
    % This function calculates the Wiener index of a given tree
    
    % Inputs:
        % T: (Matix) Adjacency matrix of the tree T
        
    % Outputs:
        % W_ind: (Integer) Wiener index of T
        % D: (Matrix) Distances between each vertex of T

    
    G = sparse(T); % converts a sparse or full matrix to sparse form
    
    D = graphallshortestpaths(G, 'directed', false); % finds all off the shortest paths in given graph
    
    M = sum(sum(D)); % calculates total distances 
    
    W_ind = M/2;

end