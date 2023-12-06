function X = shiStatExclOutlier(X,Type,Value)

% returns data after excluding outliers (assigns NaN to outliers)
% 
% X = shiStatExclOutlier(X,Type,Value)
% 
%   X     - n-by-1 or 1-by-n vector
%   Type    -   'SD'  exclude if  | value - Mean | > SD * Value (Value must be positive)
%           -   'abs' exclude if  | value - Mean | > Value (Value must be positive)
%           -   'max' exclude if  Value > max (default: )
%           -   'min' exclude if  Value < min
% 
%    ###########
% by Zhenhao Shi @ 2015-1-28
%    ###########

if isempty(X)
    return;
end

if length(size(X)) > 2 || min(size(X)) ~= 1
    error('shiStatExclOutlier only works on 1xn or nx1 matrix');
end

if nargin < 2 %default
    Type = 'Up';
    Value = 5000;
%     warning('!!!! Running an example of X < 5000 !!!!');
end

switch Type
    case {'sd','Sd','SD','Std','std','standard deviation'}
        if Value <=0
            error('SD must be above 0');
        end
        M = nanmean(X);
        SD = nanstd(X);
        Up = M+Value*SD;
        Lo = M-Value*SD;
        index = X>Up | X<Lo;
        X(index) = NaN;
    case {'absolute','Absolute','Abs','abs'}
        if Value <=0
            error('Radius must be above 0');
        end
        M = nanmean(X);
        Up = M+Value;
        Lo = M-Value;
        index = X>Up | X<Lo;
        X(index) = NaN;
    case {'up','Up','upper','Upper','Max','max','maximun','Maximun'}
        index = X>Value;
        X(index) = NaN;
    case {'lo','Lo','Low','low','lower','Lower','Min','min','minimun','Minimun'}
        index = X<Value;
        X(index) = NaN;
end