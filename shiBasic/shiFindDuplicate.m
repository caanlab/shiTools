function [idx,num] = shiFindDuplicate(MAT)

% identifies identical elements
%
% [idx,num] = shiFindDuplicate(MAT)
% idx: one value for each element in MAT; elements that that share a positive integer in idx are the same
% num: number of unique duplicates ( max(idx) )
%
%    ###########
% by Zhenhao Shi @ 2018-3-1
%    ###########

idx = zeros(size(MAT(:)));

cnt = 0;

for i = 1:length(idx)
    
    if idx(i) > 0
        continue;
    end
    
    temp_idx = shiFindDuplicate_isequal(MAT(:),MAT(i));
    
    if sum(temp_idx)>1
        cnt = cnt+1;
        idx(temp_idx) = cnt;
    end
    
end

num = max(idx);


function IX = shiFindDuplicate_isequal(ALL,ONE)

IX = false(size(ALL));

for i = 1:numel(IX)
    if isequal(ALL(i),ONE)
        IX(i) = true;
    end
end