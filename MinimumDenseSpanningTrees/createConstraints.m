function [A, b] = createConstraints(listOfEdges)
    % This function creates the constraint that the resulting vector should
    % contain at least one edge for each vertex in the graph (See Method 2 in the paper).
    
    % Input:
        % listOfEdges: (Matrix) The list of edges
        
    % Outputs:
        % A: (Matrix) Left hand side of the constraint above (Ax < b) to be used for GA
        % b: (Matrix) Right hand side of the constraint above (Ax < b) to be used for GA
    
    m = size(listOfEdges,1);
    
    numOfNodes = max(max(listOfEdges(:,1:2)));

    A = zeros(numOfNodes,m);
    b = -1*ones(numOfNodes,1);

    for i = 1:numOfNodes
        ind = find(listOfEdges(:,1) == i | listOfEdges(:,2) == i);
        A(i,ind) = -1;
    end

end