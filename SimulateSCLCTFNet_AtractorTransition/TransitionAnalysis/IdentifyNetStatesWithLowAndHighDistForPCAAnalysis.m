% Network states when distance is <= 6 and >= 40 for PCA or TSNE analysis.

clc
clear
tic

numTransit = 697;
load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

networkStates_withDistBelow6 = [];
networkStates_withDistAbove40 = [];

for i = 1:numTransit
    
    networkStates = csvread(['..\TransitionSimulations_Parallel\NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_networkStates_' num2str(i) '.csv']);
    LCSDist = transitionLCS_Hamming_Cosine_Stats{i}(:,1);
    
    
    indLowDist = find(LCSDist <= 6);
    indHighDist = find(LCSDist >= 40);
    
    tempLowDistStates = networkStates(indLowDist,:);
    tempHighDistStates = networkStates(indHighDist,:);
    
    networkStates_withDistBelow6 = [networkStates_withDistBelow6; tempLowDistStates];
    networkStates_withDistAbove40 = [networkStates_withDistAbove40; tempHighDistStates];
    
end

%csvwrite('NonNe_to_NE_networkStates_withDistanceBelow6_ad_Above40_forClustering_PCA.csv', [networkStates_withDistBelow6;networkStates_withDistAbove40])

toc
    
