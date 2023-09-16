function LCS_dist = LCS(s1, s2)

% Longest common subsequence metric as defined by CESS H. ELZINGA in
% 'Sequence analysis: metric representations of categorical time series'

    seq_difference = s1 - s2;
    matching_entity_index = find(seq_difference == 0);
    
    LCS_dist = 2 * (size(s1,2) - size(matching_entity_index,2));

end