function OutputFile = shiMoveFileToRoot(Root,FileFullPath,DestFolder,ReplaceFileSepWith)
%
% copies files under a root folder and/or its subdirectories to a single folder
%
% e.g. shiMoveFileToRoot('C:\X',{'C:\X\Y\a.nii','C:\X\Z\b.nii'},'D:','_')
%       copies the two input files in such a way:
%           C:\X\Y\a.nii  -->  D:\Y_a.nii
%           C:\X\Z\b.nii  -->  D:\Z_b.nii


if ~exist('ReplaceFileSepWith','var') || isempty(ReplaceFileSepWith)
    ReplaceFileSepWith = '_';
end

Root = char(Root);
if size(Root,1) ~= 1
    error('Root must be a single string');
end

while isequal(Root(end),filesep)
    Root = Root(1:end-1);
end
while isequal(DestFolder(end),filesep)
    DestFolder = DestFolder(1:end-1);
end
    
FileFullPath = cellstr(char(FileFullPath));

if ~isdir(Root)
    error('Root folder %s does not exist',Root);
end

File_pt = char(FileFullPath);
File_pt1 = File_pt(:,1:length(Root));
File_pt2 = File_pt(:,length(Root)+1:end);

if ~isequal(Root,char(unique(cellstr(File_pt1)))) || ~isequal(filesep,unique(File_pt2(:,1)))
    error('One or more files in FileFullPath do not have root folder Root');
end

File_pt2_new = strrep(cellstr(File_pt2(:,2:end)),filesep,ReplaceFileSepWith);

OutputFile = shiStrConcat(shiMkdir(DestFolder),filesep,File_pt2_new);

for i = 1:length(File_pt2_new)
    copyfile(FileFullPath{i},OutputFile{i});
end




    