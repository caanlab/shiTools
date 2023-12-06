function [r,t] = shiStatP2R(p,N)

% converts p value and sample size to correlation coefficient and t statistic
%
% zhenhao shi

p = 1-p/2;
t = tinv(p,N-2);
r = sqrt((t.*t)./(t.*t+(N-2)));
