% Analyze the dense spanning SCLC networks
clc
clear

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
size(nets)

nets = unique(nets,'rows');

[m, n] = size(nets);

%% Degree and Edge Analysis:

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


%% Plot the Results:
figure
[~, adj_SCLCNet] = MST(SCLCNetEdges);

G = graph(adj_SCLCNet);
G.Nodes.NodeColors = nodeDegrees(:,end);
p = plot(G);
p.NodeCData = G.Nodes.NodeColors;
colorbar
p.MarkerSize = 10;

names = {};

names{1} = 'ASCL1';names{2} = 'BCL3';names{3} = 'CEBPD';names{4} = 'EBF1';
names{5} = 'ELF3';names{6} = 'FLI1';names{7} = 'FOXA1';names{8} = 'FOXA2';
names{9} = 'GATA4';names{10} = 'GFI1B';names{11} = 'HES1';names{12} = 'ISL1';
names{13} = 'KLF2';names{14} = 'MITF';names{15} = 'MYC';names{16} = 'MYCN';
names{17} = 'NEUROD1';names{18} = 'NEUROD2';names{19} = 'NR0B1';names{20} = 'NR0B2';
names{21} = 'OLIG2';names{22} = 'POU2F3';names{23} = 'RARG';names{24} = 'RBPJ';
names{25} = 'RCOR2';names{26} = 'REST';names{27} = 'SIX5';names{28} = 'SMAD4';
names{29} = 'SOX11';names{30} = 'STAT6';names{31} = 'TCF3';names{32} = 'TCF4';
names{33} = 'TEAD4';names{34} = 'YAP1';names{35} = 'ZNF217';

p.NodeLabel = names;

for i=1:numNodes
    highlight(p,i,'MarkerSize',5 + ceil(nodeDegrees(i,end)))
end

% EdgeWeights:
% temp = tril(edgeWeights);
% [ind1, ind2] = find(temp ~= 0);
% 
% G.Edges.Weight = G.Edges.Weight*0.1;
% 
% for i = 1:size(ind1,1)
%     
%     idx = find((G.Edges.EndNodes(:,2) == ind1(i)) & (G.Edges.EndNodes(:,1) == ind2(i)));
%     
%     G.Edges.Weight(idx,1) = 1 + 2*(temp(ind1(i), ind2(i)))/max(max(temp(ind1(i), ind2(i))));
% end
% 
% 
% G.Edges.LWidths = G.Edges.Weight;
% p.LineWidth = G.Edges.LWidths;

% p.NodeColor = 'red';

% p.EdgeColor = 'k';
% p.LineWidth = 2;

%%
figure

p2 = plot(G);
p2.NodeLabel = names;
p2.NodeCData = G.Nodes.NodeColors;
colorbar
p2.MarkerSize = 10;

% EdgeWeights:
temp = tril(edgeWeights);
[ind1, ind2] = find(temp ~= 0);

G.Edges.Weight = G.Edges.Weight*0 + 0.1;

for i = 1:size(ind1,1)
    
    if(ind1(i) == 6 || ind2(i) == 6)
        idx = find((G.Edges.EndNodes(:,2) == ind1(i)) & (G.Edges.EndNodes(:,1) == ind2(i)));

        G.Edges.Weight(idx,1) = 1 + 2*(temp(ind1(i), ind2(i)))/max(max(temp(ind1(i), ind2(i))));
    end
end


G.Edges.LWidths = G.Edges.Weight;
p2.LineWidth = G.Edges.LWidths;

%%
figure
p_eg = plotOptimalTree([nets(1,:),10], SCLCNetEdges, true);
p_eg.NodeLabel = names;