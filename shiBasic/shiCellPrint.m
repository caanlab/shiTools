function CELL = shiCellPrint(varargin)

% prints 2-d cell of numbers or strings or mixture of both to screen or text file
%
% CELL = shiCellPrint(CELL)
%   print to screen, return cell matrix
% CELL = shiCellPrint(fid,CELL)
%   print to file
% CELL = shiCellPrint(...,Precision)
%   Precision - precision for number print (default='%.5f')
%
%   ###############
% Zhenhao Shi @ 2018-3-20
%   ###############

if iscell(varargin{1})
    fid = 1;
    CELL = varargin{1};
    if nargin < 2
        Precision = '%.2f';
    else
        Precision = varargin{2};
    end
else
    fid = varargin{1};
    CELL = varargin{2};
    if nargin < 3
        Precision = '%.5f';
    else
        Precision = varargin{3};
    end
end


if numel(size(CELL)) ~= 2
    error('CELL must be one- or two-D');
end


[nRow,nColumn] = size(CELL);

for i = 1:nRow
    for j = 1:nColumn
        if isnumeric(CELL{i,j})
            CELL{i,j} = num2str(CELL{i,j},Precision);
        elseif ~ischar(CELL{i,j})
            error('CELL{%d,%d} is neither numeric or character',i,j);
        end
    end
end


for i = 1:nRow
    fprintf(fid,[repmat('%s\t',1,nColumn-1),'%s\n'],CELL{i,:});
end
