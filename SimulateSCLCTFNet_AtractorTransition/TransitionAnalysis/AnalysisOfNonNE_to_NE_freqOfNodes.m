% Analysis of NonNE to NE transitions - Frequency of reactions
clc
clear

numTransit = 697;

freq = zeros(35,numTransit);

for i = 1:numTransit
    
    tractReactions = csvread(['..\TransitionSimulations_Parallel\NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_reactionTrack_' num2str(i) '.csv']);
    [m, n] = size(tractReactions);
    
    if(m>1)
        s = m;
    else
        s = n;
    end

    for j = 1:s
        
        freq(tractReactions(j),i) = freq(tractReactions(j),i) + 1;
        
    end
    
    freq(:,i) = freq(:,i)/s;
    
end

stem(sum(freq,2)/numTransit)

% csvwrite('Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_reactionWeightedFreq.csv',freq)