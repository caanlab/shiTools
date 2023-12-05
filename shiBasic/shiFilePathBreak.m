function Part = shiFilePathBreak(F)

% returns each level of directories of one or more files

F = cellstr(char(F));

if numel(F)>1
    xPart = cell(numel(F),1);
    Len = nan(numel(F),1);
    for i = 1:numel(F)
        xPart{i} = shiFilePathBreak(F{i});
        Len(i) = length(xPart{i});
    end
    Part = cell(numel(F),max(Len));
    for i = 1:numel(F)
        Part(i,1:Len(i)) = xPart{i};
    end
    return;
end

F = F{1};

[p,n,e] = fileparts(F);
Part = {[n,e]};
while ~isempty(p)
    [p,n,e] = fileparts(p);
    if isempty([n,e])
        n = p;
        p = '';
    end
    Part = [{[n,e]},Part]; %#ok<AGROW>
end


