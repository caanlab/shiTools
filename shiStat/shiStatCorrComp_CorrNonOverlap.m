function [ZPF,p] = shiStatCorrComp_CorrNonOverlap(Input1,Input2)

% compares correlated nonoverlapping correlations using Raghunathan, Rosenthal & Rubin (1996 Psychol Methods 1:178-83) procedure
%
%   It should only work for repeated-measures design, i.e. each case have
%   four values on variables X, Y, A, and B, and correlation between X and
%   Y and that between A and B are to be compared.
% 
% 
% [ZPF,p] = shiStatCorrComp_CorrNonOverlap(VariablePair_1,VariablePair_2)
%   VariablePair_1  - [X,Y]
%   VariablePair_2  - [A,B]
%   ZPF - z statistic
%   p   - significance
% 
% [ZPF,p] = shiStatCorrComp_CorrNonOverlap(CorrelationMatrix,N)
%   deals with inputs of correlation coefficients and sample size
% 
%   CorrelationMatrix - output of: corr([X,Y,A,B])
%   N - sample size
% 
%  ####################
% by Zhenhao Shi @ 2015-1-6
%  ####################


if isequal(size(Input1),size(Input2)) && size(Input1,2) == 2 && size(Input1,1) > 3
    Var = [Input1(:,1), Input2, Input1(:,2)]; % compare column 14 vs. 23 of Var, which are Var1 vs. Var2
    Var = Var(all(~isnan(Var),2),:);
    n = length(Input1);
    R = corrcoef(Var);
    r12 = R(1,2);
    r13 = R(1,3);
    r14 = R(1,4);
    r23 = R(2,3);
    r24 = R(2,4);
    r34 = R(3,4);
elseif size(Input1,1)==size(Input1,2) && sum(sum(abs(Input1-Input1')))<10^-6 && isequal(diag(Input1),ones(4,1)) && all(abs(Input1(:))<=1) && Input2>3 && Input2==int32(Input2)
    n = Input2;
    R = Input1;
    r12 = R(3,1);
    r13 = R(4,1);
    r14 = R(2,1);
    r23 = R(4,3);
    r24 = R(3,2);
    r34 = R(4,2);
else
    error('wrong input');
end


k = (r12 - r24*r14) * (r34 - r24*r23) + (r13 - r12*r23) * (r24 - r12*r14) +	(r12 - r13*r23) * (r34 - r13*r14) + (r13 - r14*r34) * (r24 - r34*r23);

% PF = (r14-r23)*sqrt(n) / sqrt( (1-r14^2)^2 + (1-r23^2)^2 - k );

z14 = atanh(r14);
z23 = atanh(r23);

ZPF = sqrt((n-3)/2) * (z14-z23) / sqrt( 1 - k/(2*(1-r14^2)*(1-r23^2)) );
p = 2*normcdf(-abs(ZPF));