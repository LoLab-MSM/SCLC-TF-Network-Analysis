function Hamm_dist = HammingDistance(s1, s2)

% Longest common subsequence metric as defined by CESS H. ELZINGA in
% 'Sequence analysis: metric representations of categorical time series'

    seq_difference = s1 - s2;
    mimmatching_entity_index = find(seq_difference ~= 0);
    
    Hamm_dist = size(mimmatching_entity_index,2);

end