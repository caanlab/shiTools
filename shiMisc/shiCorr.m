function [r] = shiCorr(x)

p = size(x,2);
r = nan(p);

% for i = 1:p-1
%     h(i,i+1:p) = corr(x(:,i),x(:,i+1:p))>thres;
% end

chunksize = 5000;
chunk = num2cell(1:chunksize:p);

for i = 1:length(chunk)
    chunk{i} = chunk{i}(1):min(p,chunk{i}(1)+chunksize-1);
end

for i = 1:length(chunk)-1
    for j = i:length(chunk)
        r(chunk{i},chunk{j}) = corr(x(:,chunk{i}),x(:,chunk{j}));
        r(chunk{j},chunk{i}) = r(chunk{i},chunk{j})';
    end
end

for i = 1:length(chunk)
    [r(chunk{i},chunk{i})] = corr(x(:,chunk{i}));
end

