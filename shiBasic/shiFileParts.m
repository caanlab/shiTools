function [Path,FileName,Ext] = shiFileParts(PathFile)

% performs "fileparts" for cellstring and returns cellstrings: paths, filenames, extentions
%
% [Path,FileName,Ext] = shiFileParts(PathFile)
%
%    ###########
% by Zhenhao Shi @ 2022-07-14
%    ###########


if isempty(PathFile)
    Path = {};
    FileName = {};
    Ext = {};
    return;
end

PathFile = cellstr(char(PathFile));

Path = cell(size(PathFile));
FileName = cell(size(PathFile));
Ext = cell(size(PathFile));

for i = 1:length(PathFile)
    [Path{i,1},FileName{i,1},Ext{i,1}] = fileparts(PathFile{i});
end