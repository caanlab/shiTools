function value = shiIf(Criterion, TrueValue, FalseValue)

% serves as a shortcut for simple if-else statements
%
% value = shiIf(Criterion, TrueValue, FalseValue)
%   returns the value as TrueValue if Criterion==TRUE, and as FalseValue if
%   Criterion==FALSE.
%   Note that both TrueValue and FalseValue must exist, otherwise use
%   traditional if-else methods.
% 
% Example 1:
%   value = shiIf(IQ<60, 'Idiot', 'NotIdiot');
%
% Example 2:
%   value = shiIf(Age<0, 0, Age);
%
%    ###########
% by Zhenhao Shi @ 2014-12-23
%    ###########

if Criterion
    value = TrueValue;
else
    value = FalseValue;
end