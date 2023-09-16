% Prepare cytoscape plotting files for some DSTs
clc
clear

networkID = 1;

SCLCNetEdges = xlsread('SCLCNetEdges.xlsx', 'C1:D239');
listOfEdges = SCLCNetEdges;    

res1 = csvread('The_dense_spanning_trees_SumSqPow560.csv');
res2 = csvread('The_dense_spanning_trees_SumSqPow576.csv');
res3 = csvread('The_dense_spanning_trees_Wiener1566.csv');
res4 = csvread('The_dense_spanning_trees_Wiener1592.csv');

nets_DST = [res1(:,1:end-1); res2(:,1:end-1); res3(:,1:end-1); res4(:,1:end-1)];
nets_DST = unique(nets_DST,'rows');

[m, n] = size(nets_DST);

numNodes = max(max(SCLCNetEdges(:,1:2)));

indOfSelectedEdges = find(nets_DST(networkID,:) == 1);

selectedEdges = [SCLCNetEdges(indOfSelectedEdges,:), ones(sum(nets_DST(networkID,:)),1)];

[adj_MST_soln, adj_G_soln] = MST(selectedEdges);

edgeStatus = [listOfEdges, zeros(size(listOfEdges,1),1)];

for i = 1:size(listOfEdges,1)
    
    sourceNode = listOfEdges(i,1);
    targetNode = listOfEdges(i,2);
        
    if(adj_MST_soln(targetNode, sourceNode) == 1)
        
        idx = find((listOfEdges(1:i-1,1) == targetNode) & (listOfEdges(1:i-1,2) == sourceNode));
        
        if(size(idx,1) == 0)
            edgeStatus(i,3) = 1;
        end
    end
        
end

filename = sprintf('AnExampleFoundDST_netID%d.xlsx',networkID);
xlswrite(filename,edgeStatus);

%% Plot of the network networkID

figure

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

p_eg = plotOptimalTree([nets_DST(networkID,:),10], SCLCNetEdges, true);
p_eg.NodeLabel = names;