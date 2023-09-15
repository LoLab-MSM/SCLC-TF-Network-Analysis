function [sumsq, D] = sumPowDegrees(T, p)

    % This function calculates the S_p of a given tree (See the paper).
    
    % Inputs:
        % T: (Matix) Adjacency matrix of the tree T
        % p: (Float) Power value for vertex degrees
        
    % Outputs:
        % W_ind: (Integer) Wiener index of T
        % D: (Matrix) Distances between each vertex of T
        
        
    [m, n] = size(T);

    G = graph(T);
    D = distances(G);

    deg = zeros(1,m);
    for i = 1:m
        counter = 0;
        for j = 1:n
            if(T(i,j) == 1)
                counter = counter + 1;
            end
        end
        deg(i) = counter;
    end

    sumsq = (sum(deg.^p)); 

end