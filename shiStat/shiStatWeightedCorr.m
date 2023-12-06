function [wR,N] = shiStatWeightedCorr(x,y,w,PositiveWeight,varargin)

% returns weighted correlation coefficient
%
% zhenhao shi

if size(w,1) == 1
    w=w';
end

if size(x,1) ~= size(y,1) || size(y,1) ~= size(w,1) || size(w,2) ~= 1
    error('unmatched input dimensions');
end


if nargin < 4 || PositiveWeight
    x = x(w>0,:);
    y = y(w>0,:);
    w = w(w>0,:);
end


pnames = {'type', 'rows'};
dflts  = {'p', 'a'};
[type,rows] = internal.stats.parseArgs(pnames,dflts,varargin{:});
typeChoices = {'pearson', 'spearman'};
rowsChoices = {'all', 'complete', 'pairwise'};
if ischar(rows)
    i = find(strncmpi(rows,rowsChoices,length(rows)));
    if isscalar(i)
        rows = rowsChoices{i}(1);
    else
        error('unkonw ''rows''');
    end
end
if ischar(type)
    i = find(strncmpi(type,typeChoices,length(type)));
    if isscalar(i)
        type = typeChoices{i}(1);
    else
        error('unkonw ''type''');
    end
end


if isequal(rows,'c')
    x = x(all(~isnan([x,y,w]),2));
    y = y(all(~isnan([x,y,w]),2));
    w = w(all(~isnan([x,y,w]),2));
end


wR = nan(size(x,2),size(y,2));
N = wR;

for i = 1:size(x,2)
    for j = 1:size(y,2)

        xx = x(:,i);
        yy = y(:,j);
        ww = w;

        if isequal(rows,'p')
            xx = xx(all(~isnan([x(:,i),y(:,j),w]),2));
            yy = yy(all(~isnan([x(:,i),y(:,j),w]),2));
            ww = ww(all(~isnan([x(:,i),y(:,j),w]),2));
        end

        if any(isnan([xx;yy;ww]))
            wR(i,j) = NaN;
            N(i,j) = 0;
            continue;
        end

        if isequal(type,'s')
            xx = tiedrank(xx);
            yy = tiedrank(yy);
        end

        wMxx = sum(xx.*ww)/sum(ww);
        wMyy = sum(yy.*ww)/sum(ww);
        wVxx = sum(ww.*(xx-wMxx).^2)/sum(ww);
        wVyy = sum(ww.*(yy-wMyy).^2)/sum(ww);
        wCov = sum(ww.*(xx-wMxx).*(yy-wMyy))/sum(ww);
        wR(i,j) = wCov/sqrt(wVxx*wVyy);
        N(i,j) = length(ww);

    end
end