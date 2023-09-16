% Identifiy increasing patterns

clc
clear
tic
load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

[m, n] = size(transitionLCS_Hamming_Cosine_Stats);

allIncreasingPatterns_NonNE2NE = {};
filteredIncreasingPatterns_NonNE2NE = {};
increasingPatternStats_NonNE2NE = {};

for t=1:n
    transition = transitionLCS_Hamming_Cosine_Stats{t};

    [m2, n2] = size(transition);

    trackIncreasePatterns = {};

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

            if(transition(i,1) >= transition(current,1))

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
        
        if((current - initStart > 8) && (size(idx,2) ~= size(tempDistValue,2)) && (max(tempDistValue) - min(tempDistValue) >= 10)) % filter the length of decreasing paths here!

            trackIncreasePatterns{counter} = [tempIndex; tempPattern; tempDistValue];

            counter = counter + 1;
        end

    end

    allIncreasingPatterns_NonNE2NE{t} = trackIncreasePatterns;
    
    % Further filtering the patters
    [m3,n3] = size(trackIncreasePatterns);

    filteredIncreasingPatterns = {};
    counter2 = 1;
    for i = 1:n3

        pattern = trackIncreasePatterns{i};

        if(pattern(3,end) >= 30)
            filteredIncreasingPatterns{counter2} = pattern;
            counter2 = counter2 + 1;
        end

    end
    
    filteredIncreasingPatterns_NonNE2NE{t}= filteredIncreasingPatterns;

    % Buraya decreasing pattern istatistiklerini cikaran kod yaz. Hagi node
    % daha cok etkilemis bu decrease'i.
    % Stats of nodes appearing in the decreasing patterns
    statsOfNodes = zeros(1,35);
    [m4,n4] = size(trackIncreasePatterns);

    for i = 1:n4

        pattern = trackIncreasePatterns{i};
        if(pattern(2,1) == -1)
            pattern = pattern(:,2:end);
        end
        statsOfNodes(1,pattern(2,:)) = statsOfNodes(1,pattern(2,:)) + 1;

    end
    
    increasingPatternStats_NonNE2NE{t} = statsOfNodes;

end

% Total Frequency of nodes in the decreasing patterns in the transition simulations 
totalFreqOfNodesInIncreasingPatterns = zeros(1,35);

for t = 1:n
    freqNodes = increasingPatternStats_NonNE2NE{t};
    
    totalFreqOfNodesInIncreasingPatterns = totalFreqOfNodesInIncreasingPatterns + freqNodes;
    
end
    
maxFreq = max(totalFreqOfNodesInIncreasingPatterns);

normalizedFreqOfNodesInDecreasingPatterns = totalFreqOfNodesInIncreasingPatterns/maxFreq;

% Save the results to .mat file

save increasingDistancePatternsAnalysisResults_NonNE2NE_2.mat filteredIncreasingPatterns_NonNE2NE increasingPatternStats_NonNE2NE totalFreqOfNodesInIncreasingPatterns normalizedFreqOfNodesInDecreasingPatterns

toc    