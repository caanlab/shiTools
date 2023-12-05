function FileName = shiFullFileName(varargin)

% returns a n-by-1 cell containing file and folder names with full paths that match wildcard identifier given by the string SuperWildCard
%
% FileName = shiFullFileName(Wildcard1, Wildcard2, ...)
% FileName = shiFullFileName(fullfile(Wildcard1, Wildcard2, ...))
% 
% Example 1: *.img files/folders in Sub* folders (full path)
%   FileName = shiFullFileName(pwd,'Sub*','*.img');
% 
% Example 2: con_0001.img files/folders in Sub* folders (full path)
%   FileName = shiFullFileName(fullfile(pwd,'Sub*','con_0001.img'));
% 
% Example 3: *PPI* files/folders (in current folder)
%   FileName = shiFullFileName('*PPI*');
%
%    ###########
% by Zhenhao Shi @ 2020-03-08
%    ###########

if nargin == 1
    SuperWildCard = varargin{1}; % single input
    if isempty(fileparts(SuperWildCard)) % check of full path
        SuperWildCard = fullfile(pwd,SuperWildCard);
    end
else
    SuperWildCard = fullfile(varargin{:}); % multiple input, concatenated
end

if iscell(SuperWildCard) % for iteration purpose
    FileName = [];
    for i = 1:length(SuperWildCard)
        FileName = [FileName;shiFullFileName(SuperWildCard{i})];
    end
    return
end

if isequal(SuperWildCard(end),filesep) % remove unnecessary ending filesep
    FileName = shiFullFileName(SuperWildCard(1:end-1));
    return;
end

Part = shiFilePathBreak(SuperWildCard);
% if length(Part{1}) > 1 && isequal(Part{1}(end-1:end),':\') % remove ending filesep in disk label in windows, e.g. change C:\ to C:
%     Part{1} = Part{1}(1:end-1);
% end
% 
FileName = shiFullFileName_xpd(Part);


function FileName = shiFullFileName_xpd(Part)
if isempty(Part)
    FileName = {};
    return;
elseif size(Part,2) == 1 % one column, job already done
    FileName = Part;
    return;
else
    if size(Part,1) > 1 % multiple columns, multiple rows, then iteratively process each row
        FileName = {};
        for i = 1:size(Part,1)
            FileName = [FileName;shiFullFileName_xpd(Part(i,:))];
        end
        return;
    end
end

% multiple columns, one row, then process first two columns
Part_1 = shiFullFileName_dir(Part{1},Part{2});
Part = [Part_1,repmat(Part(3:end),size(Part_1,1),1)];
FileName = shiFullFileName_xpd(Part); 


function p = shiFullFileName_dir(p1, p2)
px = shiFileName(p1,p2);
p = fullfile(fullfile(p1),px);