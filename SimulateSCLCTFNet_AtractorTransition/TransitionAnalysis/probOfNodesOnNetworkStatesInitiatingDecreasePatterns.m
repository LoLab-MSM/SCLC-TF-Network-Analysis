% Compute probability of nodes being ON in the network states initiating filtered decrease patterns.

clc
clear
tic

numTransit = 697;
load('decreasingDistancePatternsAnalysisResults_NonNE2NE.mat')

probNodesON = zeros(1,27);
counter = 0;

for i = 1:numTransit
    
    networkStates = csvread(['..\TransitionSimulations_Parallel\NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_networkStates_' num2str(i) '.csv']);
    decreasingPatternsPerTrans = filteredDecreasingPatterns_NonNE2NE{i}; % this is a cell for a transition
    
    [m, n] = size(decreasingPatternsPerTrans); 
    
    for j = 1:n
        
        pattern = decreasingPatternsPerTrans{j};
        
        [m2, n2] = size(pattern);
        
        counter = counter + 1;
        probNodesON = probNodesON + networkStates(pattern(1,1),:);
        
    end
    
end

probNodesON = probNodesON/counter;

% Plot the results

name = {};

for i = 1:27
    
    name{i} = TF_name(i);
    
end

f = figure;  
f.Position = [100 100 700 500]; 
stem(probNodesON, '-.h','MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0.6350 0.0780 0.1840], 'LineWidth', 1.5, 'MarkerSize',10)
ylabel('Probability', 'FontSize', 14)
xlabel('Transcription Factors', 'FontSize', 14)
title('Probability of TFs being ON in the network states initiating distance decreasing patterns', 'FontSize', 14)
xlim([0, 28])
ax = gca;
set(gca, 'XTick', 1:27)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;        

toc
