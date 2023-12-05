function Files = shiFileName_Recur(RootDir,WildCard,verbose)

% recursively searches directories and their subdirectories for files that match the WildCard
%
% Files = shiFileName_Recur(RootDir,FileFilter)
% 
%   RootDir - a string of root directory to be searched
%   FileFilter - a string for file pattern, may contain wildcard
%   Files - a cellstring of files found
%
%    ###########
% by Zhenhao Shi @ 2020-8-3
%    ###########

if ~exist('verbose','var') || isempty(verbose)
    verbose = false;
end

Files = {};

if iscell(RootDir)
    for i = 1:length(RootDir)
        Files = [Files;shiFileName_Recur(RootDir{i},WildCard,verbose)];
    end
    return;
end

x = dir(fullfile(RootDir,WildCard));
nme = {x.name}';
filt = [x.isdir]';
x = nme(~filt);
x = fullfile(RootDir,x);
if verbose
    fprintf('% 5d files found in : %s \n',length(x),RootDir);
end
Files = [Files;x];

xDir = shiFolderName(fullfile(RootDir,'*'));
if ~isempty(xDir)
    Files = [Files;shiFileName_Recur(fullfile(RootDir,xDir),WildCard,verbose)];
end
