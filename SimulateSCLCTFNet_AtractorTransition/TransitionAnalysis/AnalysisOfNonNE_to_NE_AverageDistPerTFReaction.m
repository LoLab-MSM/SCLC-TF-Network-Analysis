% Computes the average distance per TF given a transition simulation
clc
clear

load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

[m, n] = size(transitionLCS_Hamming_Cosine_Stats);

transition = transitionLCS_Hamming_Cosine_Stats{600};

[m2, n2] = size(transition);

averageDistPerTFReaction = zeros(35,1);

for i = 1:35
    
    ind = find(transition(:,4) == i);
    
    avgDist = sum(transition(ind,1))/size(ind,1);
    
    averageDistPerTFReaction(i,1) = avgDist;
    
end

%csvwrite('Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_distanceReductionAalysis.csv',transitionLCS_HammingStats)