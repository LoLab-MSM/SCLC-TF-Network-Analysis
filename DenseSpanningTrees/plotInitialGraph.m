function plotInitialGraph(adj_Init)
    % This function plots the initial given graph.
    
    % Input:
        % adj_Init: (Matrix) Adjacency matrix of the given graph

    H = graph(adj_Init);
    p = plot(H);
    p.NodeColor = 'red';
    p.MarkerSize = 10;
    p.EdgeColor = 'k';
    p.LineWidth = 2;
    
    title('The Given Graph')

end