function [vertices, cycle] = cycle_check(vertices, inserted_edge)

    % Reference : http://www.mathworks.com/matlabcentral/fileexchange/13457?kruskal?algorithm/content//MST Kruskal/iscycle.m

    % Inputs:
        % vertices: (Matrix) set of vertices in the graph
        % inserted_edge: (Row of listOfEdges), one edge we insert

    % Outputs:
        % vertices: (Matrix) The new set of vertices
        % cycle: (Integer) 1 a cycle is detected, 0 otherwise

    g = max(vertices) + 1;
    
    cycle = 0;
    
    n = length(vertices);

    if((vertices(inserted_edge(1)) == 0) && (vertices(inserted_edge(2)) == 0))
        vertices(inserted_edge(1)) = g;
        vertices(inserted_edge(2)) = g;
    elseif(vertices(inserted_edge(1)) == 0)
        vertices(inserted_edge(1)) = vertices(inserted_edge(2));
    elseif(vertices(inserted_edge(2)) == 0)
        vertices(inserted_edge(2)) = vertices(inserted_edge(1));
    elseif(vertices(inserted_edge(1)) == vertices(inserted_edge(2)))
        cycle = 1;
        return
    else
        m = max(vertices(inserted_edge(1)), vertices(inserted_edge(2)));
        for i = 1:n
            if vertices(i) == m
                vertices(i) = min(vertices(inserted_edge(1)), vertices(inserted_edge(2)));
            end
        end
    end
end