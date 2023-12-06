function [dPrime,AValue] = shiStatSdtSensitivity(Hit,FA)

% converts "hit" and "false alarm" rates to d' (d-prime) and A values
% 
% [dPrime,AValue] = shiStatSdtSensitivity(Hit,FA)
%   calculate sensitivity indices d' and A, based on:
%       1) Hit - hit rate
%       2) FA  - false alarm
% 
% AValue: Zhang, J., & Mueller, S. T. (2005). A note on ROC analysis and
%   non-parametric estimate of sensitivity. Psychometrika, 70(1), 203-212.
% 
% 
% by Zhenhao Shi @ 2015-1-6


dPrime = norminv(Hit) - norminv(FA);



AValue = nan(size(Hit));

for i = 1:numel(Hit)
    
    if Hit(i) >= FA(i)
        H=Hit(i);
        F=FA(i);
        if F<=0.5 && 0.5<=H
            A = 3/4 + (H-F)/4 - F*(1-H);
        elseif F<=H && H<0.5
            A = 3/4 + (H-F)/4 - F/(4*H);
        else
            A = 3/4 + (H-F)/4 - (1-H)/(4*(1-F));
        end;
    elseif Hit(i) < FA(i)
        A = 1/2; % ??
    else
        A = NaN;
    end
    
    AValue(i) = A;
    
end;
    
    