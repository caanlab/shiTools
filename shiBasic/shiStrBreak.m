function Part = shiStrBreak(Str,StrSep)

% breaks string(s) into parts

Str = cellstr(char(Str));

if ~exist('StrSep','var') || isempty(StrSep)
    StrSep = '_';
end

if ~ischar(StrSep)% || length(StrSep) ~= 1
    error('string separator StrSep must be a single character');
end

StringScan = ['%[^',StrSep,']%*[',StrSep,']'];

for i = 1:length(Str)
    while isequal(StrSep,Str{i}(1))
        Str{i} = Str{i}(2:end);
    end
    xPart =  textscan(Str{i},StringScan);
    xPart = xPart{1};
    Part(i,1:length(xPart)) = xPart;
end



