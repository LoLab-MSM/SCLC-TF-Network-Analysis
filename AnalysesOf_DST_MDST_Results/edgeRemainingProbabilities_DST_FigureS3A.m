clc
clear

% Edge remaining probabilities

% Stem plot the edge probabilities. X axis is the edges, y axis is the
% probabilities of them remaining in DSTs and MDSTs

% Create a new script for node names

%% DSTs %%

SCLCNetEdges = xlsread('SCLCNetEdges.xlsx', 'C1:D239');
listOfEdges = SCLCNetEdges;    
[~, adj_G_Init] = MST(listOfEdges);

res1 = csvread('The_dense_spanning_trees_SumSqPow560.csv');
res2 = csvread('The_dense_spanning_trees_SumSqPow576.csv');
res3 = csvread('The_dense_spanning_trees_Wiener1566.csv');
res4 = csvread('The_dense_spanning_trees_Wiener1592.csv');

nets_DST = [res1(:,1:end-1); res2(:,1:end-1); res3(:,1:end-1); res4(:,1:end-1)];
nets_DST = unique(nets_DST,'rows');

[m, n] = size(nets_DST);

numNodes = max(max(SCLCNetEdges(:,1:2)));
edgeProbs = zeros(numNodes);

for i = 1:m
    
    indOfSelectedEdges = find(nets_DST(i,:) == 1);

    selectedEdges = [SCLCNetEdges(indOfSelectedEdges,:), ones(sum(nets_DST(i,:)),1)];

    [adj_MST_soln, adj_G_soln] = MST(selectedEdges);
    
    ind = find(adj_MST_soln == 1);
    
    edgeProbs(ind) = edgeProbs(ind) + 1;
    
end

edgeProbs = edgeProbs/m;
edgeRemainingProbs = zeros(size(listOfEdges,1),3);
idx = 1;

for i = 1:numNodes
    for j = i+1:numNodes
        
        if(adj_G_Init(i,j) == 1)
            sourceNodeName = nodeName(j);
            targetNodeName = nodeName(i);
            
            edgeRemainingProbs(idx,1) = j;
            edgeRemainingProbs(idx,2) = i;
            edgeRemainingProbs(idx,3) = edgeProbs(i,j);
            idx = idx + 1;
        end
        
    end
end

% plot the results
numEdges2Plot = 80;

maxProb = max(edgeRemainingProbs(:,3),[],2);
[val, index] = sort(maxProb, 'descend');
sorted_edgeRemainingProbs = edgeRemainingProbs(index,:);

s = find(sorted_edgeRemainingProbs(:,1) ~= 0);
s = size(s,1);

edgeLabels = {};
for i = 1:numEdges2Plot
    sourceNodeName = nodeName(sorted_edgeRemainingProbs(i,1));
    targetNodeName = nodeName(sorted_edgeRemainingProbs(i,2));
    edgeLabels{i} = [sourceNodeName ' - ' targetNodeName];
end

figure('Position', get(0, 'Screensize'))

stem(1:numEdges2Plot,sorted_edgeRemainingProbs(1:numEdges2Plot,3), 'filled', 'LineStyle','-.','MarkerSize', 9, 'MarkerFaceColor','r');

ax = gca;
set(gca, 'XTick', 1:s)
set(gca, 'XTickLabel', edgeLabels)
ax.XTickLabelRotation = 90;
ylim([0 1.2])
ylabel('Probability of an edge remaining in the DSTs', 'FontSize', 14)

xlabel('Edges', 'FontSize', 14)

    
    