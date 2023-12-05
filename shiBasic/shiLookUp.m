function [a2b,b2a,EqualityMatrix] = shiLookUp(a,b)

% returns conversion rules between a and b such that a(a2b) == b(:), and b(b2a) == a(:)
% 
% [a2b,b2a,EqualityMatrix] = shiLookUp(a,b)
%
% in the case of non-existing elements, NaN is used
% in the case of duplicate elements, the first occurance is used
% 
%    ###########
% by Zhenhao Shi @ 2016-12-21
%    ###########


EqualityMatrix = false(numel(a),numel(b));
for ka = 1:numel(a)
    for kb = 1:numel(b)
        EqualityMatrix(ka,kb) = isequal(a(ka),b(kb));
    end
end

a2b = nan(numel(b),1);
b2a = nan(numel(a),1);

for kb = 1:numel(b)
    f = find(EqualityMatrix(:,kb));
    if isempty(f)
        fprintf('  WARNING: b(%d) = a(nothing)\n',kb);
    elseif numel(f)>1
        fprintf('  WARNING: b(%d) = a(many)     (many = %s)\n',kb,[num2str(f(:)','%d,'),sprintf('\b')]);
        a2b(kb) = f(1);
    else
        a2b(kb) = f;
    end
end

for ka = 1:numel(a)
    f = find(EqualityMatrix(ka,:));
    if isempty(f)
        fprintf('  WARNING: a(%d) = b(nothing)\n',ka);
    elseif numel(f)>1
        fprintf('  WARNING: a(%d) = b(many)     (many = %s)\n',ka,[num2str(f(:)','%d,'),sprintf('\b')]);
        b2a(ka) = f(1);
    else
        b2a(ka) = f;
    end
end


        