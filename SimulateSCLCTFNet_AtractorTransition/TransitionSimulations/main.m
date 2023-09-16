% main file to run Booleabayes rules for SCLC TF network.

    clc
    clear
    tic
%% Simulation Settings:
    numIter = 1e6;
    numTransition = 1000;

%% Read the probabity lookup tables:

    rules = fileread('rules_3_probsOnly.txt');
    rules = cellstr(rules);
    % split at newline character
    rules = cellfun(@(newline) strsplit(newline, '\n'), rules, 'UniformOutput', false);
    % reposition split values into individual cells
    % horizontal concatenation
    rules = [rules{:}];
    rules = cellfun(@str2num, rules,'un',0).';

%% Network states (Wooten et al., 2019):

    State_NE_1 = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0]; 
    State_NE_2 = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 

    State_NEv1_1 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]; 
    State_NEv1_2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0];

    State_NEv2_1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  
    State_NEv2_2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  
    State_NEv2_3 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  
    State_NEv2_4 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 

    State_NonNE_1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1];
    State_NonNE_2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1];
    
    attractors = [State_NE_1; State_NE_2; State_NEv1_1; State_NEv1_2; State_NEv2_1; ...
                  State_NEv2_2; State_NEv2_3; State_NEv2_4; State_NonNE_1; State_NonNE_2];

%% Call boolebayes function:
    
    initialState = [State_NonNE_1'; zeros(8,1)];
    
    % run booleabayes:
    for i=1:numSim
        
        [networkStates, trackReactions, transitionStatus, message] = booleabayesSCLC(numIter, initialState, rules, attractors);
        
        networkStates = networkStates';
    
        if transitionStatus == 1 % NonNE to NE transition occured
            disp(message)
            s = size(trackReactions,2);
            %h = heatmap(networkStates(1:s,1:27))
            csvwrite(['NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_networkStates_' num2str(i) '.csv'], networkStates(1:s+1,1:27))
            csvwrite(['NonNE_to_NE\Transition_from_NonNE_to_NE_ASCL1_1_FLI1_1_reactionTrack_' num2str(i) '.csv'], trackReactions)

        elseif(transitionStatus == 2)
            s = size(trackReactions,2);
            csvwrite(['NonNE_to_NEv1\Transition_from_NonNE_to_NEv1_ASCL1_1_FLI1_1_networkStates_' num2str(i) '.csv'], networkStates(1:s+1,1:27))
            csvwrite(['NonNE_to_NEv1\Transition_from_NonNE_to_NEv1_ASCL1_1_FLI1_1_reactionTrack_' num2str(i) '.csv'], trackReactions)

        elseif(transitionStatus == 3)
            s = size(trackReactions,2);
            csvwrite(['NonNE_to_NEv2\Transition_from_NonNE_to_NEv2_ASCL1_1_FLI1_1_networkStates_' num2str(i) '.csv'], networkStates(1:s+1,1:27))
            csvwrite(['NonNE_to_NEv2\Transition_from_NonNE_to_NEv2_ASCL1_1_FLI1_1_reactionTrack_' num2str(i) '.csv'], trackReactions)
        end
    end
    
    toc