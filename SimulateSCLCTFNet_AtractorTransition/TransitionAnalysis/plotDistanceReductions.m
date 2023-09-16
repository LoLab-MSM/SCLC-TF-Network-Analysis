% plot distance reductions
clc
clear

load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

[m, n] = size(transitionLCS_Hamming_Cosine_Stats);

transition1 = transitionLCS_Hamming_Cosine_Stats{1};

[m2, n2] = size(transition1);

plot(1:m2, transition1(:,1), '-.o')
ylabel('Distance to the Target State', 'FontSize', 14)
xlabel('Asynchronous Iterations ', 'FontSize', 14)