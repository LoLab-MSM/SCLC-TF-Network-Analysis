% frequency of TFs that cause a step decrease (or increase) in distance at
% each iteration

clc
clear
tic
load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

[m, n] = size(transitionLCS_Hamming_Cosine_Stats);

freqForTFsCausingDecreaseInIterations = zeros(1,35);
freqForTFsCausingIncreaseInIterations = zeros(1,35);
freqForTFsCausingNoChangeInIterations = zeros(1,35);


for t=1:n
    transition = transitionLCS_Hamming_Cosine_Stats{t};

    [m2, n2] = size(transition);

    for i = 2:m2
        
        if(transition(i,1) < transition(i-1,1))
            freqForTFsCausingDecreaseInIterations(1,transition(i,4)) = freqForTFsCausingDecreaseInIterations(1,transition(i,4)) + 1;
        elseif(transition(i,1) > transition(i-1,1))
            freqForTFsCausingIncreaseInIterations(1,transition(i,4)) = freqForTFsCausingIncreaseInIterations(1,transition(i,4)) + 1;
        elseif(transition(i,1) == transition(i-1,1))
            freqForTFsCausingNoChangeInIterations(1,transition(i,4)) = freqForTFsCausingNoChangeInIterations(1,transition(i,4)) + 1;
        end
    end
    
end

total = freqForTFsCausingDecreaseInIterations + freqForTFsCausingIncreaseInIterations + freqForTFsCausingNoChangeInIterations;
% normalized frequencies
    
maxFreqDec = max(freqForTFsCausingDecreaseInIterations);

normFreqForTFsCausingDecreaseInIterations = freqForTFsCausingDecreaseInIterations./total;


maxFreqInc = max(freqForTFsCausingIncreaseInIterations);

normFreqForTFsCausingIncreaseInIterations = freqForTFsCausingIncreaseInIterations./total;

normFreqForTFsCausingNoChangeInIterations = freqForTFsCausingNoChangeInIterations./total;


% Plot the results

name = {};

for i = 1:27
    
    name{i} = TF_name(i);
    
end

f = figure;  
f.Position = [100 100 700 500]; 
stem(freqForTFsCausingDecreaseInIterations(1,1:27), '-.h','MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0.6350 0.0780 0.1840], 'LineWidth', 1.5, 'MarkerSize',10)
ylabel('Number of Iterations', 'FontSize', 14)
xlabel('Transcription Factors', 'FontSize', 14)
title('TF Reactions Causing a Decrease in the Distance Between NonNE and NE Subtypes', 'FontSize', 14)
xlim([0, 28])
%ylim([0, 1])
ax = gca;
set(gca, 'XTick', 1:27)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;


f = figure;  
f.Position = [100 100 700 500]; 
stem(freqForTFsCausingIncreaseInIterations(1,1:27), '-.h','MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0.6350 0.0780 0.1840], 'LineWidth', 1.5, 'MarkerSize',10)
ylabel('Number of Iterations', 'FontSize', 14)
xlabel('Transcription Factors', 'FontSize', 14)
title('TF Reactions Causing an Increase in the Distance Between NonNE and NE Subtypes', 'FontSize', 14)
xlim([0, 28])
%ylim([0, 1])
ax = gca;
set(gca, 'XTick', 1:27)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;

f = figure;  
f.Position = [100 100 700 500]; 
stem(freqForTFsCausingNoChangeInIterations(1,1:27), '-.h','MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0.6350 0.0780 0.1840], 'LineWidth', 1.5, 'MarkerSize',10)
ylabel('Number of Iterations', 'FontSize', 14)
xlabel('Transcription Factors', 'FontSize', 14)
title('TF Reactions that do not Change the Distance Between NonNE and NE Subtypes', 'FontSize', 14)
xlim([0, 28])
%ylim([0, 1])
ax = gca;
set(gca, 'XTick', 1:27)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;
toc    