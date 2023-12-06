function OUT = shiTools(STR)
%
% summarizes shiTools functions


if nargin == 0
    OUT = fileparts(which('shiTools'));
    if nargout == 0
        help(fullfile(shiTools))
        help(fullfile(shiTools,'shiBasic'))
        help(fullfile(shiTools,'shiStat'))
        help(fullfile(shiTools,'shiSpm'))
        b=[
            dir(fullfile(OUT,'shiBasic','shi*'));
            dir(fullfile(OUT,'shiStat','shi*'));
            dir(fullfile(OUT,'shiSpm','shi*'));
            dir(fullfile(OUT,'shiMisc','shi*'));
            ];
        c={b.date}';
        d=sort(cellstr(datestr(c,'yyyy-mm-dd HH:MM:SS')));
        fprintf('\n shiTools: last updated at %s\n',d{end});
        return;
    end
else
    F_1 = textscan(help(fullfile(shiTools,'shiBasic')),'%[^\n\r]%*[\n\r]','Whitespace','','HeaderLines',1);
    F_2 = textscan(help(fullfile(shiTools,'shiStat')),'%[^\n\r]%*[\n\r]','Whitespace','','HeaderLines',1);
    F_3 = textscan(help(fullfile(shiTools,'shiSpm')),'%[^\n\r]%*[\n\r]','Whitespace','','HeaderLines',1);
    ALL = [F_1{1};F_2{1};F_3{1}];
    OUT = char(ALL(shiStrIncl(ALL,STR,false)));
end
    
    