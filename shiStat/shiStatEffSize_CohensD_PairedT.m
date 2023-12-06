function d = shiStatEffSize_CohensD_PairedT(Dat1,Dat2)

% computes effect size as Cohen's d for paired t-tests
% 
%   Dat1   - n-by-p matrix, where n is the number of cases, p
%            is the number of variables
%   Dat2   - n-by-p matrix
%    d     - 1-by-p vector of Cohen's d values for each variable
% 
% 
% by Zhenhao Shi @ 2020-4-7
% 
% 

if ~isequal(size(Dat1),size(Dat2))
    error('unmatched number of cases/variables');
end

d = shiStatEffSize_CohensD_OneT(Dat1-Dat2);


