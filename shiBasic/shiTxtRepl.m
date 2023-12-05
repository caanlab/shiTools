function shiTxtRepl(OldFile,Old,New,NewFile,isRegExp)

% performs "replace text1 with text 2" and "save as a new file" for a text file
%
% shiTxtRepl(OldFile,Old,New,NewFile)
%   replaces Old strings in OldFile by New strings and save as NewFile.
%   OldFile is a string of the file name of the to-be-replaced text file.
%   Old is a string or a cell array of regular expressions for strings to
%   be replaced. New is a string or a cell array of strings with the same
%   size as Old, replacing correponding elements of Old. Replaced content
%   is saved in the file with name NewFile.
% 
% Example:
%   shiTxtRepl('Sub01_batch.m','Sub01','Sub02','Sub02_batch.m')
%
%    ###########
% by Zhenhao Shi @ 2019-8-4
%    ###########

if ischar(Old)
    Old = {Old};
end
if ischar(New)
    New = {New};
end

if ~exist('isRegExp','var') || isempty(isRegExp)
    isRegExp = true;
end

if isRegExp
    REPL_FUN = @regexprep;
else
    REPL_FUN = @strrep;
end

[~,o1,o2] = fileparts(OldFile);
[~,n1,n2] = fileparts(NewFile);
fprintf('replacing: %s -> %s\n',[o1,o2],[n1,n2]);

if numel(Old) ~= numel(New)
    error('unmatched size');
end

Text = cellstr(readlines(OldFile));
% OLD VERSION FOR READING
% fid = fopen(OldFile);
% Text = textscan(fid,'%[^\n\r]%*[\n\r]','Whitespace','');
% fclose(fid);
% Text = Text{1};

for i = 1:length(Old)
    Text = REPL_FUN(Text,Old{i},New{i});
end
Text = regexprep(Text,'[\n\r]','');


fid = fopen(NewFile,'wt+');
for i = 1:length(Text)
    fprintf(fid,'%s\n',Text{i});
end
fclose(fid);