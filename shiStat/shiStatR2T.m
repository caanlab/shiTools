function [t,p] = shiStatR2T(r,N)

% converts correlation coefficient to t-statistic
%
% zhenhao shi

t = r.*sqrt((N-2)./(1-r.^2));
p = tcdf(t,N-2);
p = min(p,1-p)*2;
