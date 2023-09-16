% Statistical Analysis of the DSTs and MDSTs
clc
clear

%% Initial graph:
% node degrees
% node degree distributions

SCLCNetEdges = xlsread('SCLCNetEdges.xlsx', 'C1:D239');

listOfEdges = SCLCNetEdges;
[m, n] = size(listOfEdges);
    
[~, adj_G_Init] = MST(listOfEdges);
%plotInitialGraph(adj_G_Init)

numNodes = max(max(SCLCNetEdges(:,1:2)));
nodeDegrees = zeros(numNodes, 1);
nodeDegrees(:,1) = sum(adj_G_Init,2);    

figure
subplot(1,3,1)
h1 = histogram(nodeDegrees, 15);
ylim([0,15])
xlim([0,30])
xlabel('Node degree', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)
title('The initial network node degree distribution', 'FontSize', 14)
h1.FaceColor = [0 0.4470 0.7410];


%% DSTs-MSTs:
% average node degrees
% node degree distributions

%%%% DST %%%%:
SCLCNetEdges = xlsread('SCLCNetEdges.xlsx', 'C1:D239');

res1 = csvread('The_dense_spanning_trees_SumSqPow560.csv');
res2 = csvread('The_dense_spanning_trees_SumSqPow576.csv');
res3 = csvread('The_dense_spanning_trees_Wiener1566.csv');
res4 = csvread('The_dense_spanning_trees_Wiener1592.csv');

nets_DST = [res1(:,1:end-1); res2(:,1:end-1); res3(:,1:end-1); res4(:,1:end-1)];
nets_DST = unique(nets_DST,'rows');

[m, n] = size(nets_DST);

numNodes = max(max(SCLCNetEdges(:,1:2)));
nodeDegrees = zeros(numNodes, m + 1);
edgeWeights = zeros(numNodes, numNodes);

for i = 1:m
    
    indOfSelectedEdges = find(nets_DST(i,:) == 1);

    selectedEdges = [SCLCNetEdges(indOfSelectedEdges,:), ones(sum(nets_DST(i,:)),1)];

    [adj_MST_soln, adj_G_soln] = MST(selectedEdges);
    
    nodeDegrees(:,i) = sum(adj_MST_soln,2);
    
    ind = find(adj_MST_soln == 1);
    
    edgeWeights(ind) = edgeWeights(ind) + 1;
    
end
        
nodeDegrees(:,m+1) = sum(nodeDegrees(:,1:m),2)/m;

subplot(1,3,2)
h2 = histogram(nodeDegrees(:,end), 15);
ylim([0,35])
xlim([0,30])
xlabel('Node degree', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)
title('DST average node degree distribution', 'FontSize', 14)
h2.FaceColor = [0.9290 0.6940 0.1250];

%%%% MDST %%%%:
SCLCNetEdges = csvread('SCLCnetwork_edgeSourceTarget_withProbs.csv');
SCLCNetEdges(:,3) = abs(SCLCNetEdges(:,3));
SCLCNetEdges(:,3) = 1 - SCLCNetEdges(:,3);
idx = find(SCLCNetEdges(:,3) ~= 1);
SCLCNetEdges = SCLCNetEdges(idx,:);

res1 = csvread('The_minimum_dense_spanning_trees_SumSq1_530.csv');
res2 = csvread('The_minimum_dense_spanning_trees_SumSq2_530.csv');
res3 = csvread('The_minimum_dense_spanning_trees_SumSq3_496.csv');
res4 = csvread('The_minimum_dense_spanning_trees_Wiener1_1724.csv');
res5 = csvread('The_minimum_dense_spanning_trees_Wiener2_1724.csv');
res6 = csvread('The_minimum_dense_spanning_trees_Wiener3_1724.csv');
res7 = csvread('The_minimum_dense_spanning_trees_Wiener4_1724.csv');


nets = [res1(1:end,1:end-1); res2(:,1:end-1); res3(:,1:end-1); res4(:,1:end-1); res5(:,1:end-1); res6(:,1:end-1); res6(:,1:end-1)];
nets = unique(nets,'rows');

[m, n] = size(nets);

numNodes = max(max(SCLCNetEdges(:,1:2)));
nodeDegrees = zeros(numNodes, m + 1);
edgeWeights = zeros(numNodes, numNodes);

for i = 1:m
    
    indOfSelectedEdges = find(nets(i,:) == 1);

    selectedEdges = SCLCNetEdges(indOfSelectedEdges,:);

    [adj_MST_soln, adj_G_soln] = MST(selectedEdges);
    
    nodeDegrees(:,i) = sum(adj_MST_soln,2);
    
    ind = find(adj_MST_soln == 1);
    
    edgeWeights(ind) = edgeWeights(ind) + 1;
    
end
        
nodeDegrees(:,m+1) = sum(nodeDegrees(:,1:m),2)/m;

subplot(1,3,3)
h3 = histogram(nodeDegrees(:,end), 15);
ylim([0,35])
xlim([0,30])
xlabel('Node degree', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)
title('MDST average node degree distribution', 'FontSize', 14)
h3.FaceColor = [0.4940 0.1840 0.5560];