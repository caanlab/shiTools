function DirMade = shiMkdir(varargin)

% makes new directory(s), if not existed
%
% DirMade = shiMkdir(Dir)
% 
%   Dir - string or cell array of strings for directory(s) to be created
%   DirMade - full path of new directory(s)
% 
%    ###########
% by Zhenhao Shi @ 2018-7-17
%    ###########

if nargin == 1
    Dir = varargin{1};
else
    Dir = fullfile(varargin{:});
end

isCharDir = ischar(Dir);

Dir = cellstr(char(Dir));
DirMade = cell(size(Dir));

for i = 1:length(Dir)
    try
        if ~isdir(char(shiFullFileName(Dir{i})))
            mkdir(Dir{i});
        end
        DirMade{i} = char(shiFullFileName(Dir{i}));
    catch
        error('failed to create: %s\n',Dir{i});
        DirMade{i} = '';
    end
end

if isCharDir
    DirMade = char(DirMade);
end