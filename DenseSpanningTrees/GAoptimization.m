function [uniqueSols2, fval] = GAoptimization(listOfEdges, wantDense, isObjFunctMinimized, populationSize, useParallel)
    % This function runs Matlab's GA based on user specified settings.
    
    % Inputs:
        % listOfEdges: (Matrix) The list of edges where each row represent an edge with its weight.
        % wantDense: (Boolean) true if DST is searched, false if SST is searched.
        % isObjFunctMinimized: (Boolean) true if the objective function is minimized for the desired tree, false otherwise.
        % populationSize: (Integer) The initial population size.
        % useParallel: (Boolean) true if parallel computing tool is used.
        
    % Outputs:
        % uniqueSols: (Matrix) The matrix where each row is a vector h representing the found tree.
        % fval: (Integer or Float) The best objective function value that GA could find.
    
    
    listEdges = listOfEdges;
    
    m = size(listEdges,1);
    
    % Creating the constraints by calling createConstraints() function
    [A, b] = createConstraints(listEdges);
    
    % Run 1:
    if(1e4 < populationSize)
        options = optimoptions('ga', 'PopulationSize', 1e4, 'UseParallel', useParallel);
    else
        options = optimoptions('ga', 'PopulationSize', populationSize, 'UseParallel', useParallel);
    end
    
    
    [x,fval,exitflag,output,population,scores] = ga(@objFunct,m,A,b,[],[],zeros(1,m),ones(1,m),[],1:1:m,options);
    
    solutions1 = [population, scores];
    ind = find(solutions1(:,end) == fval);
    bestSolutions1 = solutions1(ind,:);
    uniqueSols1 = unique(bestSolutions1,'rows');
    
    % Run 2:
    if(size(uniqueSols1,1) < populationSize)
        options = optimoptions('ga', 'PopulationSize', populationSize, 'UseParallel', useParallel, 'InitialPopulationMatrix', uniqueSols1(:,1:end-1));
    else
        options = optimoptions('ga', 'PopulationSize', size(uniqueSols1,1), 'UseParallel', useParallel, 'InitialPopulationMatrix', uniqueSols1(:,1:end-1));
    end
    
    [x,fval,exitflag,output,population,scores] = ga(@objFunct,m,A,b,[],[],zeros(1,m),ones(1,m),[],1:1:m,options);
    
    solutions2 = [population, scores];
    
    % Some found solutions may not have the same objective value. Hence, we filter the ones with worse objective value.
    
    ind = find(solutions2(:,end) == fval);
    bestSolutions2 = solutions2(ind,:);
    uniqueSols2 = unique(bestSolutions2,'rows');

    function objValue = objFunct(selectedEdges)

        indOfSelectedEdges = find(selectedEdges == 1);

        G_selected = listEdges(indOfSelectedEdges,:);

        [adj_MST_selected, adj_G_selected] = MST(G_selected);
        
        % The following "objective function", i.e., wiener_index(), method can be replaced by another desired objective function. 
        % The code will not require extra adjustments. Just use the same input/output format like given objective functions here.
        
        [f_val, D] = wiener_index(adj_MST_selected);
        % [f_val, D] = sumPowDegrees(adj_MST_selected, 2);
        
        [i, j] = find(D == Inf); % checks whether the tree is connected or not
        if(size(i,1) > 0) % If the found tree is not connected, then we assign a very high objective value so that it won't survive in the next generation of GA
            objValue = 10e10;
        else
            % The following conditional statement is necessary because by default GA minimizes the its @objFunct. So, if, in case, the
            % maximization wants to be done, then the objective value is multiplied by -1 as (max f(x)) == (min -f(x))
            if(wantDense && isObjFunctMinimized)
                objValue = f_val;
            elseif(wantDense && (~isObjFunctMinimized))
                objValue = -f_val;
            elseif((~ wantDense) && isObjFunctMinimized)
                objValue = f_val;
            elseif((~ wantDense) && (~isObjFunctMinimized))
                objValue = -f_val;    
            end
        end
         
    end
    
end