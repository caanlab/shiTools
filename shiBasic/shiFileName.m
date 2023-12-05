function FileName = shiFileName(varargin)

% returns a n-by-1 cell containing file and/or folder names that match wildcard identifier given by the string WildCard
%
% FileName = shiFileName(WildCard)
% 
% Example 1: *.img files in current directory
%   FileName = shiFileName('*.img');
% 
% Example 2: *.img files in a given directory (absolute path)
%   FileName = shiFileName('C:\Study\ImageData\*.img');
% 
% Example 3: swra* files in a given directory (relative path)
%   FileName = shiFileName('ImageData/swra*');
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
    FileName = shiFileName(WildCard(1:end-1));
    return;
end

if ~isempty(regexp(WildCard,'*','once'))
    FileName = dir(WildCard);
    FileName = {FileName.name}';
else
    if isfolder(WildCard)
        [~,FileName1,FileName2] = fileparts(WildCard);
        FileName = {[FileName1,FileName2]};
    else
        FileName = dir(WildCard);
        FileName = {FileName.name}';
    end
end

if isempty(FileName)
    FileName = {};
else
    FileName = FileName(~shiStrIncl(FileName,{'^\.$','^\.\.$'}));
end

