% This is the main file that controls every step of finding dense and sparse spanning trees from a given graph using the proposed Method 2.
% Theoretical information on Method 2 can be found in the following paper:
    % M. Ozen, G. Lesaja, and H. Wang, Globally optimal dense and sparse spanning trees, and their applications, Statistics, Optimization & Information Computing, 8(2), pp. 328-345, 2020.

% Instructions:
    % 1. Load your list of edges. Your file should contain 2 (unweighted edges) or 3 (weighted edges) columns. Each row should contain
        % end-vertex labels of an edge. For instance, a row [1, 5, 0.5] --> means an edge that is connected between vertex 1 and vertex 5 with weight 0.5.

    % 2. Specify the objective function. In the code, Weiner Index and Sum of power(p) of degrees are included. Any desired objective function can be implemented and embedded in this code.
        % To change the objective function, go to line 46 in GAoptimization.m file.

    % 3. Specify objective function settings:
        % wantDense = true, if you want to search for dense spanning tree
        % wantDense = false, if you want to search for sparse spanning tree
        % isObjFunctMinimized = true, if the desired tree is found by minimizing the objective function
        % isObjFunctMinimized = false, if the desired tree is found by maximizing the objective function
        % Note: The dense(sparse) spanning trees can minimize(maximize) some objective functions whihe they may maximize(minimize) some other objective functions. Therefore, type of optimization should be specified.
        % DSTs minimize Weiner index whereas SSTs maximize Weiner index. The opposite is true for Sum of power(p) of degrees for p > 1.

    % 4. Genetic algorithm settings:
        % Set the populationSize to a number. A high population size increases the chance of finding better dense/sparse spanning trees while it increases the time complexity. 
        % Also, for large graphs, it is suggested to use a high population size (eg. 5e5, etc.)

        % Set useParallel = true if you want to run the code in parallel mode. Running in parallel mode will speed up the code as all cores of the computer are used.
        % Note: You should have the Matlab's Parallel Computing Toolbox installed in your system to run in parallel.

% The testcase(by default):
    % The used dataset is a randomly generated graph containing both densest(star) and sparsest (path) trees. 
    % (See Ozen, M. et al.,Globally optimal dense and sparse spanning trees, and their applications, Statistics, Optimization & Information Computing, 8(2), pp. 328-345, 2020)

    % The objective function: Wiener Index (minimize for DSTs and maximize for SSTs)

    % wantDense = true
    % isObjFunctMinimized = true

    % populationSize = 1e4; 
    % useParallel = false;
clc
clear
tic
%% -------------------------- PART I: Preparation -------------------------
% -------------------------------------------------------------------------

% 1. Load dataset
    SCLCNetEdges = xlsread('SCLCNetEdges.xlsx', 'C1:D239');
    
    listOfEdges = SCLCNetEdges;
    [m, n] = size(listOfEdges);
    
% 2. See instructions above.

% 3. Objective function settings:
    wantDense = true;
    isObjFunctMinimized = true; 

% 4. GA settings:
    populationSize = 1e5; 
    useParallel = true;

% The following step checks if the loaded dataset contains weighted edges. If it is not weighted, then we add equal weight 1 to all edges.    
    if(n == 2)
        listOfEdges = [listOfEdges, ones(m,1)];
    end
 
% The following step checks if the vertex labels start from 0. If yes, then we add 1 to each vertex label so that the vertex labels become 1 to N.    
    if(min(min(listOfEdges(:,1:2))) == 0)
        listOfEdges(:,1:2) = listOfEdges(:,1:2) + 1;
    end
       
% Visualize the given initial graph
    [~, adj_G_Init] = MST(listOfEdges);
    plotInitialGraph(adj_G_Init)

%% ------------------------ PART II: GA Optimization ----------------------
% -------------------------------------------------------------------------

    [uniqueSols, finalObjValue] = GAoptimization(listOfEdges, wantDense, isObjFunctMinimized, populationSize, useParallel);

%% ------------------------ PART III: Display and Save --------------------
% -------------------------------------------------------------------------
    figure
    plotOptimalTree(uniqueSols, listOfEdges, wantDense)
    
    finalObjValue = abs(finalObjValue)

    if(wantDense)
        csvwrite('The_dense_spanning_trees.csv', uniqueSols)
    else
        csvwrite('The_sparse_spanning_trees.csv', uniqueSols)
    end
toc