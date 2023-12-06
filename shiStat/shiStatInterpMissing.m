function Y_New = shiStatInterpMissing(Y,InterpOrder)

% replaces missing values in 1-D time series (equal interval) by interpolation
% 
% Y_New = shiStatInterpMissing(Y)
% Y_New = shiStatInterpMissing(Y,InterpOrder)
%   Y           - raw 1-D time series
%   InterpOrder - default = 'linear';
%   Y_New       - time series with missing values replaced
% 
% by Zhenhao Shi
% 


if ~isvector(Y)
    error('Y must be a vector');
end


if nargin < 2
    InterpOrder = 1;
end

if isscalar(InterpOrder)
    switch InterpOrder
        case 0
            Method = 'nearest';
        case 1
            Method = 'linear';
        case 2
            Method = 'spline';
        otherwise
            error('unknown interpolation method');
    end
else
    Method = InterpOrder;
end
                
Y_New = nan(size(Y));
Y_New(:) = interp1(find(~isnan(Y)),Y(~isnan(Y)),1:length(Y),Method);