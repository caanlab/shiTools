function d = shiStatEffSize_CohensD_OneT(Dat)

% computes effect size as Cohen's d for one-sample t-tests
% 
%   Dat   - n-by-p matrix, where n is the number of cases, p
%            is the number of variables
%   d     - 1-by-p vector of Cohen's d values for each variable
% 
% 
% by Zhenhao Shi @ 2020-4-7
% 
% 

M = nanmean(Dat);
S = nanstd(Dat);
d = M./S;


