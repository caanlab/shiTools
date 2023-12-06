function [Z,p] = shiStatCorrComp_CorrOverlap(Input1,Input2)

% compares correlated overlapping correlations using Meng, Rosenthal & Rubin (1992 Psychol Bull 111:172-5) procedure
%
%   It should only work for repeated-measures design, i.e. each case have
%   three values on variables X1, X2, and Y, and correlation between X1 and
%   Y and that between X2 and Y are to be compared.
% 
% [Z,p] = shiStatCorrComp_CorrOverlap(Var_X12,Var_Y)
%   Var_X12         - [X1,X2]
%   Var_Y           -  Y
%   Z        - z statistic
%   p        - significance
% 
% [Z,p] = shiStatCorrComp_CorrOverlap(CorrelationMatrix,N)
%   deals with inputs of correlation coefficients [r1 r2 rx] or r matrix and sample size
% 
%   CorrelationMatrix - output of: corr([X1,X2,Y])
%   N - sample size
% 
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########

if isequal(size(Input1,1),size(Input2,1)) && size(Input1,2) == 2 && size(Input1,1) > 3 && size(Input2,2) == 1
    VarX1 = Input1(:,1);
    VarX2 = Input1(:,2);
    VarY = Input2;
    ind = ~isnan(VarX1) & ~isnan(VarX2) & ~isnan(VarY);
    VarX1 = VarX1(ind);
    VarX2 = VarX2(ind);
    VarY = VarY(ind);
    n = length(VarX1);
    r1 = corr(VarX1,VarY);
    r2 = corr(VarX2,VarY);
    rx = corr(VarX1,VarX2);
elseif size(Input1,1)==size(Input1,2) && sum(sum(abs(Input1-Input1')))<10^-6 && isequal(diag(Input1),ones(3,1)) && all(abs(Input1(:))<=1) && Input2>3 && Input2==int32(Input2)
    n = Input2;
    r1 = Input1(3,1);
    r2 = Input1(3,2);
    rx = Input1(2,1);
elseif numel(Input1)==3 && Input2>3 && Input2==int32(Input2)
    n = Input2;
    r1 = Input1(1);
    r2 = Input1(2);
    rx = Input1(3);
else
    error('wrong input');
end




meanR2 = (r1*r1+r2*r2)/2;

f = (1-rx)/(2*(1-meanR2));

h = 1+meanR2*(1-f)/(1-meanR2);

z_r1 = atanh(r1);
z_r2 = atanh(r2);

Z = (z_r1-z_r2)*sqrt((n-3)/(2*h*(1-rx)));
p = 2*normcdf(-abs(Z));