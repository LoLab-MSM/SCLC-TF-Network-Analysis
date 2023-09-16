function [T, adj_G] = MST(listOfEdges)

    % This function finds the initial minimum spanning tree from a given
    % graph using Kruskal's Algorithm (See Method 2 in the paper for details).
    
    % Input:
        % listOfEdges = (Matrix) The list of edges including the weights of edges.
    
    % Outputs:
        % T: (Matrix) Adjacency matrix of the minimum spanning tree
        % adj_G: (Matrix) Adjacency matrix of the given graph


    m = size(listOfEdges,1);
    
    numNodes = max(max(listOfEdges(:,1:2)));
    adj_G = zeros(numNodes);

    % Create adjacency matrix of the given graph
    for i = 1:m
       adj_G(listOfEdges(i,1), listOfEdges(i,2)) = 1;
       adj_G(listOfEdges(i,2), listOfEdges(i,1)) = 1;
    end

    % Sort edges based on their weights by ascending order
    listOfEdges = sort(listOfEdges, 3, 'ascend');
    vertices = zeros (1, numNodes);
    T = zeros(numNodes);

    % Check if the graph has cycle after adding an edge
    for i = 1:m
       
       inserted_edge = listOfEdges(i,1:2);
       [vertices, cycle] = cycle_check(vertices, inserted_edge);
       
       if(cycle == 1)
           listOfEdges(i,:) = 0;
       end
       
    end

    % Create minimum spanning tree's adjacency matrix
    for i = 1:m
       
       if(listOfEdges(i,1:2) ~= 0)
           T(listOfEdges(i,1), listOfEdges(i,2)) = 1;
           T(listOfEdges(i,2), listOfEdges(i,1)) = 1;
       end
       
    end
end