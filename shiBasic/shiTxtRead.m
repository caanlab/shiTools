function T = shiTxtRead(Filename)

% reads the content of a text file into a cellstr
% 
% T = shiTxtRead(Filename)
%   returns the content of the text file with name Filename to a n-by-1
%   cell array of strings T, where n is the number of lines in Filename.
%
%    ###########
% by Zhenhao Shi @ 2016-7-26
%    ###########


T = cellstr(readlines(Filename));




%% OLD VERSION
% 
% fid = fopen(Filename);
% if fid<0
%     error('%s does not exist',Filename);
% end
% T = textscan(fid,'%[^\n\r]%*[\n\r]','Whitespace','');
% fclose(fid);
% 
% T = T{1};