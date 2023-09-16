% Compute probability of nodes being ON when a decrease occurs.

clc
clear
tic

numTransit = 697;
load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

probNodesON = zeros(1,27);
counter = 0;

for i = 1:numTransit
    
    networkStates = csvread(['..\TransitionSimulations_Parallel\NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_networkStates_' num2str(i) '.csv']);
    transition = transitionLCS_Hamming_Cosine_Stats{i};
    
    [m, n] = size(networkStates);
    
    prev = transition(1,1);
    for j = 2:m
        
        if(transition(j,1) < prev)
            counter = counter + 1;
            probNodesON = probNodesON + networkStates(j,:);
            prev = transition(j,1);
        else
            prev = transition(j,1);
        end
        
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
title('Probability of TFs being ON when decrease in distance occurs', 'FontSize', 14)
xlim([0, 28])
ax = gca;
set(gca, 'XTick', 1:27)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;        

toc
