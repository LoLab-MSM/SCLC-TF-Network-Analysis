function p = plotOptimalTree(uniqueSol, listOfEdges, wantDense)

    % This function plots one of the found optimal trees.
    
    % Inputs:
        % uniqueSol: (Matrix) The list of found solution vectors h
        % listOfEdges: (Matrix) The list of edges
        % wantDense: (Boolean) true if DST is searched, false if SST is searched.
    
    indOfSolnEdges = find(uniqueSol(1,1:end-1) == 1);

    G_soln = listOfEdges(indOfSolnEdges,:);

    [adj_MST_soln, adj_G_soln] = MST(G_soln);

    H = graph(adj_MST_soln);
    p = plot(H);
    p.NodeColor = 'red';
    p.MarkerSize = 10;
    p.EdgeColor = 'k';
    p.LineWidth = 2;
    
    if(wantDense)
        title('The Dense Spanning Tree')
    else
        title('The Sparse Spanning Tree')
    end
end