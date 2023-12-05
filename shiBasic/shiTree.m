function Tree = shiTree(Path)

% returns all sub-dir under a root dir
% 
% Tree = shiTree(Path)
%   returns all paths under the path given by the string Path to the n-by-1
%   cell array Tree, where n is the number of paths within Path. Tree
%   includes Path itself. If Path is a full/relative path, so are the
%   elements of Tree.
%
%    ###########
% by Zhenhao Shi @ 2015-2-27
%    ###########

Path = shiFullFileName(Path);
if numel(Path) ~= 1
    error('invalid root path');
end
Path = Path{1};

Tree = genpath(Path);
Tree = regexprep(Tree,pathsep,''';''');
eval(['Tree = {''',Tree,'''};']);
Tree = Tree(1:end-1);