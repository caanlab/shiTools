function Strucel = shiStrRepl(Strucel,Old,New,isRegExp)

% recursively searches for and replaces old strings with new strings in a cell or struct variable
%
% Strucel = shiStrRepl(Strucel,Old,New)
% Strucel = shiStrRepl(Strucel,Old,New,isRegExp)
%   Strucel - a cell or a struct variable that may contain string(s)
%   Old - old string pattern to be looked for and replaced
%   New - new string that the old one is to be replaced with
%   isRegExp - whether "Old" is a regular expression (default = true). if true, @regexprep is used. if false, @strrep is used
%
% Example: 
%   SPM_new = shiStrRepl(SPM,'C:\AnalysisFolder','/Volume/HardDrive/Analysis',false)
%   SPM_new = shiStrRepl(SPM,'\','/',false)
%   (will update directory and filesep in SPM)
%
% Zhenhao Shi 2020/1/20

if ~exist('isRegExp','var') || isempty(isRegExp)
    isRegExp = true;
end

if iscell(Strucel)
    for i = 1:numel(Strucel)
        Strucel{i} = shiStrRepl(Strucel{i},Old,New,isRegExp);
    end
    return
elseif isstruct(Strucel)
    xField = fieldnames(Strucel);
    for i = 1:numel(Strucel)
        for f = 1:numel(xField)
            Strucel(i).(xField{f}) = shiStrRepl(Strucel(i).(xField{f}),Old,New,isRegExp);
        end
    end
    return
end

if isRegExp
    REPL_FUN = @regexprep;
else
    REPL_FUN = @strrep;
end

if ischar(Strucel) || isstring(Strucel)
    Strucel = REPL_FUN(Strucel,Old,New);
end
    
