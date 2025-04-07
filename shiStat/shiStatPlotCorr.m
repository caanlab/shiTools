function [rho,H,fDots,fLines] = shiStatPlotCorr(Y,X,Group,YName,XName,Title,CorrType,verbose)
 
% gives Y-X scatterplot as well as linear trendline by groups, and prints correlation coefficient and p-values at command window
%
% shiStatPlotCorr(Y,X)
% 
%   Y    - n-by-1 vector, where n is the number of observations
%   X    - n-by-1 vector
% shiStatPlotCorr(Y,X,Group)
%   divides observations into groups based on Group, and plot and calculate
%   Y-X correlation individually.
% 
%   Group - n-by-1 vector or cell array of strings. Observations from the
%           same group have the same value in Group, where as those from a
%           different group have a different value.
%           default = ones(n,1)
% 
%   Example:
%           Y = rand(10,1);
%           X = rand(10,1);
%           Group1 = [1 1 2 2 1 2 1 1 2 1]';
%           Group2 = {'male'
%                     'female'
%                     'male'
%                     'female'
%                     'male'
%                     'male'
%                     'female'
%                     'female'
%                     'male'
%                     'male'};
%           shiStatPlotCorr(Y,X,Group1);
%           shiStatPlotCorr(Y,X,Group2);
% 
% shiStatPlotCorr(Y,X,Group,YName)
%   assign y-axis the name YName (default = 'Y');
% 
% shiStatPlotCorr(Y,X,Group,YName,XName)
%   assign x-axis the name XName (default = 'X');
% 
% shiStatPlotCorr(Y,X,Group,YName,XName,Title)
%   assign plot the name Title (default: YName - XName);
% 
% shiStatPlotCorr(Y,X,Group,YName,XName,Title,CorrType)
%   CorrType can be 'pearson' (default), 'spearman', or 'kendall'
% 
% H = shiStatPlotCorr(...)
%   returns handle of the scatter plot
% 
%    ###########
% by Zhenhao Shi @ 2024-8-1
%    ###########


if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

if ~exist('CorrType','var') || isempty(CorrType)
    CorrType = 'p';
end

Color = {'r','b','g','m','c','k'};
% Marker = {'.','*','s','^','d','x'};
Marker = {'.','.','.','.','.','.'};

Y = reshape(Y,numel(Y),1);
X = reshape(X,numel(X),1);

Y(isinf(Y)) = NaN;
X(isinf(X)) = NaN;

if ~isequal(size(Y),size(X))
    error('unmatched number of observations');
end

if nargin < 3
    Group = ones(numel(Y),1);
elseif isempty(Group)
    Group = ones(numel(Y),1);
else
    Group = reshape(Group,numel(Group),1);
end

if isnumeric(Group) || islogical(Group)
    Group = mat2cell(num2str(Group),ones(1,size(num2str(Group),1)),size(num2str(Group),2));
end

NOTNAN = intersect(find(~isnan(X)),find(~isnan(Y)));

Y = Y(NOTNAN);
X = X(NOTNAN);
Group = Group(NOTNAN);

GroupLabel = unique(Group);
GroupIndex = cell(length(GroupLabel),1);

for i = 1:length(GroupLabel)
    GroupIndex{i} = strcmp(GroupLabel{i},Group);
    if sum(GroupIndex{i})<3
        error('too few observations in Group: %s',GroupLabel{i});
    end
end

if length(GroupLabel) == 1
    Color = {'k'};
    Marker = {'.'};
elseif length(GroupLabel) <= 6
    Color = Color(1:length(GroupLabel));
    Marker = Marker(1:length(GroupLabel));
else
    Color = mat2cell(hsv(length(GroupLabel)),ones(1,length(GroupLabel)),3);
    Marker = repmat({'.'},1,length(GroupLabel));
end


rho = nan(length(GroupLabel),1);
sig = nan(length(GroupLabel),1);
P = cell(length(GroupLabel),1);
R = cell(length(GroupLabel),1);

fDots = cell(length(GroupLabel),1);
fLines = cell(length(GroupLabel),1);

for i = 1:length(GroupLabel)
    [rho(i),sig(i)] = corr(X(GroupIndex{i}),Y(GroupIndex{i}),'type',CorrType);
    N(i) = length(X(GroupIndex{i}));
    P{i} = polyfit(X(GroupIndex{i}),Y(GroupIndex{i}),1);
    R{i} = P{i}(1) .* X(GroupIndex{i}) + P{i}(2);
    fDots{i} = plot(X(GroupIndex{i}), Y(GroupIndex{i}), Marker{i}, 'MarkerEdgeColor', Color{i},'MarkerSize',20);
    hold on;
end

for i = 1:length(GroupLabel)
    if rho(i) > 0
        fLines{i} = plot([min(X(GroupIndex{i})),max(X(GroupIndex{i}))], [min(R{i}),max(R{i})], '-','Color',Color{i}, 'LineWidth', 1);
    else
        fLines{i} = plot([min(X(GroupIndex{i})),max(X(GroupIndex{i}))], [max(R{i}),min(R{i})], '-','Color',Color{i}, 'LineWidth', 1);
    end
    hold on;
end

hold off;

axis([min(X)-(max(X)-min(X))/20 max(X)+(max(X)-min(X))/20 ...
    min(Y)-(max(Y)-min(Y))/20 max(Y)+(max(Y)-min(Y))/20]);

if length(GroupLabel)>1
    set(legend([GroupLabel;GroupLabel]),'Interpreter','none','Location','best');
end

if nargin < 4
    YName = 'Y';
end
if nargin < 5
    XName = 'X';
end
if nargin < 6
    Title = [YName,' - ',XName];
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
  





