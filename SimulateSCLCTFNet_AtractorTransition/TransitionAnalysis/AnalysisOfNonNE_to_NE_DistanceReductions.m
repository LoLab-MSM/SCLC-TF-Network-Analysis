% Computes distance reductions during transition simulations based on LCS,
% Hamming, and Cosine distances. 
% !!! NO need to RUN as the results of this code saved in:
% !!! transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat!!!

clc
clear

numTransit = 697;

%State_NonNE_1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1];

transitionLCS_Hamming_Cosine_Stats = {};

for i = 1:numTransit
    
    tractReactions_temp = csvread(['..\TransitionSimulations_Parallel\NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_reactionTrack_' num2str(i) '.csv']);
    networkStates = csvread(['..\TransitionSimulations_Parallel\NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_networkStates_' num2str(i) '.csv']);
    
    [m, n] = size(networkStates);
    
    converged_NE_state = networkStates(end,:);
    
    LCS_temp = zeros(m,1);
    Hamming_temp = zeros(m,1);
    Cosine_temp = zeros(m,1);
    
    for j = 1:m
        
        LCS_temp(j,1)= LCS(converged_NE_state,networkStates(j,:));
        Hamming_temp(j,1) = HammingDistance(converged_NE_state,networkStates(j,:));
        Cosine_temp(j,1) = CosineDist(converged_NE_state,networkStates(j,:));
        
    end
        
        
    [m2, n2] = size(tractReactions_temp);
    if(m2 == 1)
        tractReactions = [-1, tractReactions_temp]';
    else
        tractReactions = [-1; tractReactions_temp];
    end
    
    transitionLCS_Hamming_Cosine_Stats{i} = [LCS_temp, Hamming_temp, Cosine_temp, tractReactions];
    
end

save transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat transitionLCS_Hamming_Cosine_Stats
%csvwrite('Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_distanceReductionAalysis.csv',transitionLCS_HammingStats)