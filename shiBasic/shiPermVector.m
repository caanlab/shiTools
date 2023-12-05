function Y = shiPermVector(V,ChunkSize)

% randomly shuffles chunks (instead of single elements) of a vector
%
% Y = shiPermVector(V,ChunkSize)
% 
%    ###########
% by Zhenhao Shi @ 2018-2-13
%    ###########


if ~isvector(V)
    error('input must be a vector');
end

if nargin < 2
    ChunkSize = 1;
end

nWin = ceil(length(V)/ChunkSize);

Index = cell(1,nWin);
Index_reduce = randi(nWin);

dd = 0;
for i = 1:nWin
    Index{i} = (dd+1):(dd+ChunkSize);
    if i == Index_reduce
        Index{i} = Index{i}(1:(ChunkSize-(ChunkSize*nWin-length(V))));
    end
    dd = Index{i}(end);
end

Index = cell2mat(Index(randperm(nWin)));
Y = V(Index);