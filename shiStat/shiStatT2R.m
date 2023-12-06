function [r,p] = shiStatT2R(t,df)

% converts correlation coefficient to t-statistic
%
% zhenhao shi

r = t./sqrt(t.^2+df);
p = tcdf(t,df);
p = min(p,1-p)*2;
