% Identifiy TFs initiating decreasing patterns

clc
clear
tic
load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

[m, n] = size(transitionLCS_Hamming_Cosine_Stats);

allDecreasingPatterns_NonNE2NE = {};
filteredDecreasingPatterns_NonNE2NE = {};
TFsInitiatingDecreasingPatternStats_NonNE2NE = {};

for t=1:n
    transition = transitionLCS_Hamming_Cosine_Stats{t};

    [m2, n2] = size(transition);

    trackDecreasePatterns = {};

    start = 1;
    counter = 1;
    while(start < m2)

        initStart = start;

        tempIndex = [];
        tempPattern = [];
        tempDistValue = [];
        tempPattern = [tempPattern, transition(start,4)];
        tempDistValue = [tempDistValue, transition(start,1)];
        tempIndex = [tempIndex, start];

        distValueAtStart = transition(start,1);

        current = start;
        for i = start+1:m2 % Check for infinite loop!

            if(transition(i,1) <= transition(current,1))

                tempIndex = [tempIndex, i];
                tempPattern = [tempPattern, transition(i,4)];
                tempDistValue = [tempDistValue, transition(i,1)];
                current = i;

            else
                break;
            end

        end

        start = current + 1;
        
        val = tempDistValue(1,1);
        idx = find(tempDistValue == val);
        
        if((current - initStart > 8) && (size(idx,2) ~= size(tempDistValue,2))) % filter the length of decreasing paths here!

            trackDecreasePatterns{counter} = [tempIndex; tempPattern; tempDistValue];

            counter = counter + 1;
        end

    end

    allDecreasingPatterns_NonNE2NE{t} = trackDecreasePatterns;
    
    % Further filtering the patters
    [m3,n3] = size(trackDecreasePatterns);

    filteredDecreasingPatterns = {};
    counter2 = 1;
    for i = 1:n3

        pattern = trackDecreasePatterns{i};

        if(pattern(3,end) < 10)
            filteredDecreasingPatterns{counter2} = pattern;
            counter2 = counter2 + 1;
        end

    end
    
    filteredDecreasingPatterns_NonNE2NE{t}= filteredDecreasingPatterns;

    % Stats of nodes appearing in the decreasing patterns
    statsOfNodesInitiatingDecrease = zeros(1,35);
    [m4,n4] = size(trackDecreasePatterns);

    for i = 1:n4

        pattern = trackDecreasePatterns{i};
        if(pattern(2,1) == -1)
            pattern = pattern(:,2:end);
        end
        
        ind = 1;
        temp = pattern(2,1);
        while(ind < size(pattern(2,:),2))
            
            if(pattern(3,ind+1) < pattern(3,ind))
                temp = pattern(2,ind+1);
                break;
            else
                ind = ind + 1;
            end
        end
            
        statsOfNodesInitiatingDecrease(1,temp) = statsOfNodesInitiatingDecrease(1,temp) + 1;

    end
    
    TFsInitiatingDecreasingPatternStats_NonNE2NE{t} = statsOfNodesInitiatingDecrease;

end

% Total Frequency of nodes in the decreasing patterns in the transition simulations 
totalFreqOfNodesInitiatingDecreasingPatterns = zeros(1,35);

for t = 1:n
    freqNodes = TFsInitiatingDecreasingPatternStats_NonNE2NE{t};
    
    totalFreqOfNodesInitiatingDecreasingPatterns = totalFreqOfNodesInitiatingDecreasingPatterns + freqNodes;
    
end
    
maxFreq = max(totalFreqOfNodesInitiatingDecreasingPatterns);

normalizedFreqOfNodesInitiatingDecreasingPatterns = totalFreqOfNodesInitiatingDecreasingPatterns/maxFreq;

% Plot the results
name = {};

for i = 1:27
    
    name{i} = TF_name(i);
    
end

f = figure;  
f.Position = [100 100 700 500]; 
stem(normalizedFreqOfNodesInitiatingDecreasingPatterns(1,1:27), '-.h','MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0.6350 0.0780 0.1840], 'LineWidth', 1.5, 'MarkerSize',10)
ylabel('Normalized Frequency', 'FontSize', 14)
xlabel('Transcription Factors', 'FontSize', 14)
title('Frequency of TF Reactions Initiating Distance Decreasing Patterns', 'FontSize', 14)
xlim([0, 27])
ax = gca;
set(gca, 'XTick', 1:28)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;

toc    