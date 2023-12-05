function Incl = shiStrIncl(Target,Badge,CaseSensitive)

% (obsolete; consider "contains" "endsWith" "startsWith") uses regular expression to determine if each element of a cellstr matches any of the given expressions
%
% Incl = shiStrIncl(Target,Badge)
%   compares each element of Target, which is a string or a cell array of 
%   strings, to all elements of Badge, which is a string or a cell array of
%   regular expression strings. If the Target element matches any of the
%   Badge elements, TRUE is returned to the correponding element of Incl, 
%   which is a logical matrix of the same size as Target.
%
% Incl = shiStrIncl(Target,Badge,CaseSensitive)
%   assign CaseSensitive as 0 to ignore case (default = 1)
% 
% Example:
%   Incl = shiStrIncl({'Sub01';'Sub02';'Tub01';'Tun02'},{'^Sub','ub01$'});
%   Incl =
%           [TRUE; TRUE; TRUE; FALSE];
%
%    ###########
% by Zhenhao Shi @ 2018-6-19
%    ###########

if nargin<3
    CaseSensitive = 1;
end

if ischar(Target)
    Target = {Target};
end

if iscell(Target)
    for i = 1:size(Target,2)
        if ~iscellstr(Target(:,i))
            Target(:,i) = cellstr(char(Target(:,i)));
        end
    end
end
if isempty(Badge) || isequal(Badge,{''}) || isequal(Badge,{[]})
    Incl = false(size(Target));
    return;
end
if ~iscellstr(Badge)
    Badge = cellstr(char(Badge));
end

Incl = zeros(size(Target))>0;

if CaseSensitive
    for i = 1:numel(Target)
        Incl(i) = 0;
        for j = 1:numel(Badge)
            if ~isempty(regexp(Target{i},Badge{j}, 'once'))
                Incl(i) = true;
                break;
            end
        end
    end
else
    for i = 1:numel(Target)
        Incl(i) = 0;
        for j = 1:numel(Badge)
            if ~isempty(regexpi(Target{i},Badge{j}, 'once'))
                Incl(i) = true;
                break;
            end
        end
    end
end    