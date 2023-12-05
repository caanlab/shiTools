function FolderName = shiFolderName(varargin)

% returns a n-by-1 cell containing folder names that match wildcard identifier given by the string WildCard.
% 
% FolderName = shiFolderName(WildCard)
%
% Example 1: Sub* folders in current directory
%   FolderName = shiFolderName('Sub*');
% 
% Example 2: Sub* folders in a given directory (absolute path)
%   FolderName = shiFolderName('C:\Study\ImageData\Sub*');
% 
% Example 3: *ppi folders in a given directory (relative path)
%   FolderName = shiFolderName('ImageData/*ppi');
%
%    ###########
% by Zhenhao Shi @ 2018-7-17
%    ###########

if nargin == 1
    WildCard = varargin{1};
else
    WildCard = fullfile(varargin{:});
end

if isequal(WildCard(end),filesep)
    FolderName = shiFolderName(WildCard(1:end-1));
    return;
end

if ~isempty(regexp(WildCard,'*','once'))
    FolderName = dir(WildCard);
    FolderName = {FolderName(cell2mat({FolderName.isdir})).name}';
%     FolderName = struct2cell(dir(WildCard));
%     FolderNameAll = FolderName(1,:)';
%     FolderIndex = cell2mat(FolderName(4,:));
%     FolderName = FolderNameAll(FolderIndex);
else
    if isdir(WildCard)
        [~,FolderName1,FolderName2] = fileparts(WildCard);
        FolderName = {[FolderName1,FolderName2]};
    else
        FolderName = {};
    end
end;

if ~isempty(FolderName)
    FolderName = FolderName(~shiStrIncl(FolderName,{'^\.$','^\.\.$'}));
end