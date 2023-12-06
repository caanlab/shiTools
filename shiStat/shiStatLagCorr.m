function [r,p,Lag,Z,N] = shiStatLagCorr(nLag,x,varargin)

% returns lagged correlation coefficients
%
% [r,p,Lag,Z,N] = shiStatLagCorr(nLag,x,varargin)
% 
%   nLag        - a non-negative integer, for maximum number of lags
%   x,varagin   - see corr.m
%   r           - 3-D lagged correlation coefficient matrix, where the
%                 third diminsion covers the range of lags (see output Lag)
%                 the matrix will also be plotted
%   p           - 3-D p value matrix
%   Lag         - -nLag:nLag, corresponding to the third dimension of other
%                 ouputs. Negative lag is x predicting y, and vice versa
%   Z           - Fisher's Z transformation of r
%   N           - 3-D matrix of sample size
%
% 
% by Zhenhao Shi @ 2015-6-21
% 


if nargin < 2
    error('Too few input');
end

if ~isscalar(nLag) || nLag < 0
    error('nLag must be a non-negative integer');
end

if (nargin < 3) || ischar(varargin{1})
    y = x;
else
    y = varargin{1};
    varargin = varargin(2:end);
end



pnames = {'type', 'rows', 'tail'};
dflts  = {'p', 'a', 'both'};
[~,rows,~] = internal.stats.parseArgs(pnames,dflts,varargin{:});
rowsChoices = {'all', 'complete', 'pairwise'};
if ischar(rows)
    i = find(strncmpi(rows,rowsChoices,length(rows)));
    if isscalar(i)
        rows = rowsChoices{i}(1);
    else
        error('unkonw ''rows'' type');
    end
end



Lag = -nLag:nLag;

r = nan(size(x,2),size(y,2),length(Lag));
p = nan(size(x,2),size(y,2),length(Lag));
N = nan(size(x,2),size(y,2),length(Lag));

for lag = 1:length(Lag)

    xx = x(max(1,1-Lag(lag)):min(end,end-Lag(lag)),:);
    yy = y(max(1,1+Lag(lag)):min(end,end+Lag(lag)),:);

    [r(:,:,lag),p(:,:,lag)] = corr(xx,yy,varargin{:});
    N(:,:,lag) = funcshi_CorrCountN(xx,yy,rows);

end

Z = atanh(r);



k = 0;
if isvector(x) && isvector(y)
    plot(Lag,squeeze(r(1,1,:)));
else
    for p1 = 1:size(x,2)
        for p2 = 1:size(y,2)
            k = k+1;
            subplot(size(x,2),size(y,2),k);
            plot(Lag,squeeze(r(p1,p2,:)));
        end
    end
end




function N = funcshi_CorrCountN(xx,yy,rows)

N = nan(size(xx,2),size(yy,2));

switch rows

    case 'a'
        N = ones(size(xx,2),size(yy,2)) * size(xx,1);

    case 'c'
        nn = sum(~isnan(sum(xx,2)+sum(yy,2)));
        N = ones(size(xx,2),size(yy,2)) * nn;

    case 'p'
        for p1 = 1:size(xx,2)
            for p2 = 1:size(yy,2)
                N(p1,p2) = sum(~isnan(xx(:,p1)+yy(:,p2)));
            end
        end

end
