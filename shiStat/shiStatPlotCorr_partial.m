function [H,xY,xX,Group] = shiStatPlotCorr_partial(Y,X,C,Group,YName,XName,Title,CorrType,verbose)

% makes group-wise partial correlation plot, similar to shiStatPlotCorr
% 
% H = shiStatPlotCorr_partial(Y,X,C,Group,YName,XName,Title,CorrType)
%   where C is the covariate(s)
%
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########


if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

if ~exist('CorrType','var') || isempty(CorrType)
    CorrType = 'Pearson';
end

Color = {'r','b','g','m','c','k'};
Marker = {'.','*','s','^','d','x'};

Y = reshape(Y,numel(Y),1);
X = reshape(X,numel(X),1);

Y(isinf(Y)) = NaN;
X(isinf(X)) = NaN;
C(isinf(C)) = NaN;

if ~isequal(size(Y),size(X))
    error('unmatched number of observations');
end

if nargin < 4
    Group = ones(numel(Y),1);
elseif isempty(Group)
    Group = ones(numel(Y),1);
else
    Group = reshape(Group,numel(Group),1);
end

if isnumeric(Group) || islogical(Group)
    Group = mat2cell(num2str(Group),ones(1,size(num2str(Group),1)),size(num2str(Group),2));
end

tmp = sum([X,Y,C],2);
NOTNAN = find(~isnan(tmp));

Y = Y(NOTNAN);
X = X(NOTNAN);
C = C(NOTNAN,:);
Group = Group(NOTNAN);

GroupLabel = unique(Group);
GroupIndex = cell(length(GroupLabel),1);

for i = 1:length(GroupLabel)
    GroupIndex{i} = strcmp(GroupLabel{i},Group);
    if sum(GroupIndex{i})<3
        error('too few observations in Group: %s',GroupLabel{i});
    end
end

if length(GroupLabel) <= 6
    Color = Color(1:length(GroupLabel));
    Marker = Marker(1:length(GroupLabel));
else
    Color = repmat({'k'},1,length(GroupLabel));
    Marker = repmat({'o'},1,length(GroupLabel));
end


rho = nan(length(GroupLabel),1);
sig = nan(length(GroupLabel),1);
P = cell(length(GroupLabel),1);
R = cell(length(GroupLabel),1);

xX = shiStatResid(X,C,true,true);
xY = shiStatResid(Y,C,true,true);

for i = 1:length(GroupLabel)
    [rho(i),sig(i)] = partialcorr(X(GroupIndex{i}),Y(GroupIndex{i}),C(GroupIndex{i},:),'Type',CorrType);
    N(i) = length(X(GroupIndex{i}));
    P{i} = polyfit(xX(GroupIndex{i}),xY(GroupIndex{i}),1);
    R{i} = P{i}(1) .* xX(GroupIndex{i}) + P{i}(2);
    plot(xX(GroupIndex{i}), xY(GroupIndex{i}), Marker{i}, 'MarkerEdgeColor', Color{i},'MarkerSize',20);
    hold on;
end

for i = 1:length(GroupLabel)
    if rho(i) > 0
        plot([min(xX(GroupIndex{i})),max(xX(GroupIndex{i}))], [min(R{i}),max(R{i})], ['-',Color{i}], 'LineWidth', 1);
    else
        plot([min(xX(GroupIndex{i})),max(xX(GroupIndex{i}))], [max(R{i}),min(R{i})], ['-',Color{i}], 'LineWidth', 1);
    end
    hold on;
end

hold off;

axis([min(X)-(max(X)-min(X))/20 max(X)+(max(X)-min(X))/20 ...
    min(Y)-(max(Y)-min(Y))/20 max(Y)+(max(Y)-min(Y))/20]);

if length(GroupLabel)>1
    set(legend([GroupLabel;GroupLabel]),'Interpreter','none');
end

if nargin < 5
    YName = 'Y';
end
if nargin < 6
    XName = 'X';
end
if nargin < 7
    Title = [YName,' - ',XName,' (partial)'];
end
set(xlabel(XName),'Interpreter','none');%,'FontName','Calibri','FontSize',14);
set(ylabel(YName),'Interpreter','none');%,'FontName','Calibri','FontSize',14,'Rotation',90);
set(title(Title),'Interpreter','none');

if verbose
    fprintf('%s\n',Title);
    for i = 1:length(GroupLabel)
        fprintf('  Group = %s : r=%+.2f, N=%d, p=%g\n', GroupLabel{i}, rho(i), N(i), sig(i));
    end
    fprintf('\n');
end

if nargout > 0
    H = gcf;
end
  





