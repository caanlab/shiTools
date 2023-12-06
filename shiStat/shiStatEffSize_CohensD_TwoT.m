function d = shiStatEffSize_CohensD_TwoT(Dat1,Dat2)

% computes effect size as Cohen's d for two-sample t-tests
% 
%   Dat1   - n1-by-p matrix, where n1 is the number of cases in group 1, p
%            is the number of variables
%   Dat2   - n2-by-p matrix, where n2 is the number of cases in group 2
%    d     - 1-by-p vector of Cohen's d values for each variable
% 
% 
% by Zhenhao Shi @ 2020-4-7
% 
% 

if size(Dat1,2)~=size(Dat2,2)
    error('unmatched number of variables');
end

M1 = nanmean(Dat1);
M2 = nanmean(Dat2);
S1 = nanstd(Dat1);
S2 = nanstd(Dat2);
N1 = sum(~isnan(Dat1));
N2 = sum(~isnan(Dat2));
S = sqrt( ((N1-1).*S1.^2+(N2-1).*S2.^2) ./ (N1+N2-2) );
d = (M1-M2)./S;


