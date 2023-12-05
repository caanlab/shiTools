function CellStr = shiStrConcat(varargin)

% connects cellstrs, strings, and numbers with singleton expansion and leading zero
%
% CellStr = shiStrConcat(varargin)
%   note:   1. trailing spaces in cellstr varagin{i} will be removed
%           2. leading zero will apply to number varagin{i} to ensure equal
%           width
%           3. equal width for cellstr varagin{i} will not be guarenteed
%
% Example:
% s = shiStrConcat({'abc';'rst';'xyz'},'_',9:11,'_',1)
% 
% s =
% 
%   3×1 cell array
% 
%     {'abc_09_1'}
%     {'rst_10_1'}
%     {'xyz_11_1'}
%
%    ###########
% by Zhenhao Shi @ 2020-7-1
%    ###########

ISEMPTY = true;
for i = 1:nargin
    if ~isempty(varargin{i})
        ISEMPTY = false;
        break;
    end
end
if ISEMPTY
    CellStr = {''};
    return;
end

ISEMPTY2 = false;
for i = 1:nargin
    if iscell(varargin{i}) && isempty(varargin{i})
        ISEMPTY2 = true;
        break;
    end
end
if ISEMPTY2
    CellStr = {};
    return;
end

SizeArgin = nan(nargin,1);
for i = 1:nargin
    if isempty(varargin{i})
        SizeArgin(i) = 1;
        continue;
    end
    if isnumeric(varargin{i})
        varargin{i} = varargin{i}(:);
        if isequal(varargin{i},int32(varargin{i})) && max(varargin{i})>0
            varargin{i} = num2str(varargin{i},['%0',num2str( max(floor(log10(max(varargin{i}))),0) +1),'d']);
        else
            varargin{i} = num2str(varargin{i}(:),['%0',num2str(size(num2str(varargin{i}(:),'%0g'),2)),'g']);
        end
    end
    if isequal(varargin{i},' ')
        varargin{i} = {' '};
    end
    if ischar(varargin{i})
        varargin{i} = cellstr(char(varargin{i}));
    end
    if iscellstr(varargin{i})
        varargin{i} = varargin{i}(:);
    end
    if ~iscellstr(varargin{i}) || numel(size(varargin{i}))>2 || size(varargin{i},2)~=1
        error('wrong input');
    end
    SizeArgin(i) = numel(varargin{i});
end

SizeArgin = unique(SizeArgin);
if numel(SizeArgin)>2
    error('incompatible length of strings');
elseif numel(SizeArgin)==2 && min(SizeArgin)~=1
    error('incompatible length of strings');
end

nRow = max(SizeArgin);
for i = 1:nargin
    if length(varargin{i}) == 1
        varargin{i} = repmat(varargin{i},[nRow,1]);
    end
end

allStr = horzcat(varargin{:});
CellStr = cell(nRow,1);
for r = 1:nRow
    CellStr{r} = horzcat(allStr{r,:});
end


