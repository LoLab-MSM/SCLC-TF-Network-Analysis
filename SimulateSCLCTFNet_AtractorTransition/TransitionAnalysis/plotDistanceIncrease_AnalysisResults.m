% plot distance reductions

clc
clear

load('transitionLCS_Hamming_Cosine_Stats_NonNE2NE.mat')

transitionLCS_HammingStats = transitionLCS_Hamming_Cosine_Stats;

[m, n] = size(transitionLCS_HammingStats);

transition1 = transitionLCS_HammingStats{1};

[m2, n2] = size(transition1);

%% Whole transition iterations
% s = stairs(1:m2, transition1(:,1), '-.o', 'LineWidth', 1);
% ylabel('LCS-based Distance to the Target State', 'FontSize', 14)
% xlabel('Asynchronous Iterations ', 'FontSize', 14)
% title('Change in LCS-based Distance in NonNE to NE Asynchronous Transition')
% set(gca, 'box', 'off');

%Alternative
x = 1:m2;
y = transition1(:,1)';
z = y;
% patch([x nan],[y nan],[z nan],[z nan], 'edgecolor', 'interp');
patch([x nan],[y nan],[z nan],[z nan], 'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
colormap summer
colorbar('XTickLabel',{'', ''},'XTick', [0, 46])
caxis ([0 46])
ylabel('LCS-based Distance to the Target State', 'FontSize', 14)
xlabel('Asynchronous Iterations ', 'FontSize', 14)
title('Change in LCS-based Distance in NonNE to NE Asynchronous Transition')
set(gca, 'box', 'off');
h=text(1,1,'NonNE Subtype');
set(h,'Rotation',-90);
h=text(0,0,'NE Subtype');
set(h,'Rotation',-90);


%% Increase pattern 1
%patternInd = 1536:1:1592;
patternInd = 67113:1:67184;


name = {};

for i = 1:size(patternInd,2)
    
    name{i} = TF_name(transition1(patternInd(1,i),4));
    
end

% figure
% stairs(patternInd, transition1(patternInd',1), '-.o', 'LineWidth', 1.5, 'MarkerSize',10)
% ylabel('LCs-based Distance to the Target State', 'FontSize', 14)
% xlabel('Asynchronous Iterations ', 'FontSize', 14)
% xlim([1469, 1537])
% ax = gca;
% set(gca, 'XTick', patternInd)
% set(gca, 'XTickLabel', name)
% set(gca, 'box', 'off');
% ax.XTickLabelRotation = -90;
% grid on

%Alternative
figure
x = patternInd;
y = transition1(patternInd',1)';
z = y;
% patch([x nan],[y nan],[z nan],[z nan], 'edgecolor', 'interp');
patch([x nan],[y nan],[z nan],[z nan], 'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
colormap summer
caxis ([0 46])
ylabel('LCs-based Distance to the Target State', 'FontSize', 14)
xlabel('Asynchronous Iterations ', 'FontSize', 14)
xlim([67112, 67185])
ylim([0, 45])
ax = gca;
set(gca, 'XTick', patternInd)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;
grid on

%% Increase pattern 2
patternInd = 298442:1:298543;

name = {};

for i = 1:size(patternInd,2)
    
    name{i} = TF_name(transition1(patternInd(1,i),4));
    
end

% figure
% stairs(patternInd, transition1(patternInd',1), '-.o', 'LineWidth', 1.5, 'MarkerSize',10)
% ylabel('LCS-based Distance to the Target State', 'FontSize', 14)
% xlabel('Asynchronous Iterations ', 'FontSize', 14)
% xlim([m2-52, m2+1])
% ax = gca;
% set(gca, 'XTick', patternInd)
% set(gca, 'XTickLabel', name)
% set(gca, 'box', 'off');
% ax.XTickLabelRotation = -90;
% grid on

%Alternative
figure
x = patternInd;
y = transition1(patternInd',1)';
z = y;
% patch([x nan],[y nan],[z nan],[z nan], 'edgecolor', 'interp');
patch([x nan],[y nan],[z nan],[z nan], 'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
colormap summer
caxis ([0 46])
ylabel('LCS-based Distance to the Target State', 'FontSize', 14)
xlabel('Asynchronous Iterations ', 'FontSize', 14)
xlim([298441, 298544])
ax = gca;
set(gca, 'XTick', patternInd)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;
grid on

%% Frequency of nodes in increasing patterns

load('increasingDistancePatternsAnalysisResults_NonNE2NE.mat')

name = {};

for i = 1:27
    
    name{i} = TF_name(i);
    
end

f = figure;  
f.Position = [100 100 700 500]; 
stem(totalFreqOfNodesInIncreasingPatterns(1,1:27), '-.h','MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0.6350 0.0780 0.1840], 'LineWidth', 1.5, 'MarkerSize',10)
ylabel('Number of Occurance', 'FontSize', 14)
xlabel('Transcription Factors', 'FontSize', 14)
title('Occurance of TF Reactions in Distance Increasing Patterns', 'FontSize', 14)
xlim([0, 28])
ax = gca;
set(gca, 'XTick', 1:27)
set(gca, 'XTickLabel', name)
set(gca, 'box', 'off');
ax.XTickLabelRotation = -90;