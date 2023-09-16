function cosine_dist = CosineDist(s1, s2)

% Longest common subsequence metric as defined by CESS H. ELZINGA in
% 'Sequence analysis: metric representations of categorical time series'
    
    temp = [s1;s2];
    cosine_dist = pdist(temp,'cosine');

end