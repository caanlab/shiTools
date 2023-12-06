function h = shiCorrThres(x,varargin)
%
% returns logical matrix indicating whether correlation coefficients are above certain a threshold
%
% can deal with a large number of variables
%
% shiCorrThres(x,thres)
% shiCorrThres(x,thres,'type',CorrType,'rows',RowType,'tail',TailType) [see corr.m]
% shiCorrThres(x,y,thres)
% shiCorrThres(x,y,thres,'type',CorrType,'rows',RowType,'tail',TailType) [see corr.m]
%
% input: thres - cells of 2-vectors specifying range, [L,H), i.e. true if L <= r < H
%                e.g.: thres = {[-0.5, 0.5]; [0.5, Inf]; [-Inf, 0.5]};
%
% Zhenhao Shi 2020-5-10


chunk_size_factor = 7;

if isnumeric(varargin{1}) && size(varargin{1},1)==size(x,1)
    y = varargin{1};
    thres = varargin{2};
    varargin = varargin(3:end);
    chunksize = round(sqrt(size(x,2))*chunk_size_factor);
    chunksizey = round(sqrt(size(y,2))*chunk_size_factor);
    corrXX = false;
    px = size(x,2);
    py = size(y,2);
else
    thres = varargin{1};
    varargin = varargin(2:end);
    chunksize = round(sqrt(size(x,2))*chunk_size_factor);
    corrXX = true;
    px = size(x,2);
    py = [];
end

if ~iscell(thres)
    error('thres must be cells of 2-vectors specifying range: [L,H).\n   e.g. thres = {[-0.5, 0.5]; [0.5, Inf]; [-Inf, 0.5]};');
end

if ~isempty(varargin) && isequal(varargin{end},'thresOK')
    varargin = varargin(1:end-1);
else
    MaxNumThres = shiChkMem(px,py,numel(thres));
    fprintf('max number of thresholds worked on at the same time: %d\n', MaxNumThres);
    fprintf('number of chunks for x: %d\n', chunksize);
    if ~corrXX
        fprintf('number of chunks for y: %d\n', chunksizey);
    end
end

if exist('MaxNumThres','var') && MaxNumThres < numel(thres)
    h = cell(numel(thres),1);
    for it = 1:ceil(numel(thres)/MaxNumThres)
        indBeg = (it-1) * MaxNumThres+1;
        indEnd = min( it * MaxNumThres, numel(thres));
        if ~corrXX
            varargin = [{y},{thres(indBeg:indEnd)},varargin(:)'];
        else
            varargin = [{thres(indBeg:indEnd)},varargin(:)'];
        end
        h(indBeg:indEnd) = shiCorrThres(x,varargin{:},'thresOK');
    end
    return
end


chunkx = num2cell(1:chunksize:px);
for i = 1:length(chunkx)
    chunkx{i} = chunkx{i}(1):min(px,chunkx{i}(1)+chunksize-1);
end

if ~corrXX

    chunky = num2cell(1:chunksizey:py);
    for j = 1:length(chunky)
        chunky{j} = chunky{j}(1):min(py,chunky{j}(1)+chunksizey-1);
    end

    h = cell(numel(thres),1);
    for k = 1:numel(thres)
        h{k} = true(px,py);
    end

    for i = 1:length(chunkx)
        for j = 1:length(chunky)
            coef = corr(x(:,chunkx{i}),y(:,chunky{j}),varargin{:});
            for k = 1:numel(thres)
                h{k}(chunkx{i},chunky{j}) = coef >= thres{k}(1) & coef < thres{k}(2);
            end
        end
    end

else

    h = cell(numel(thres),1);
    for k = 1:numel(thres)
        h{k} = true(px,px);
    end

    for i = 1:length(chunkx)-1
        for j = i:length(chunkx)
            coef = corr(x(:,chunkx{i}),x(:,chunkx{j}),varargin{:});
            for k = 1:numel(thres)
                h{k}(chunkx{i},chunkx{j}) = coef >= thres{k}(1) & coef < thres{k}(2);
                h{k}(chunkx{j},chunkx{i}) = h{k}(chunkx{i},chunkx{j})';
            end
        end
    end
    for i = 1:length(chunkx)
        coef = corr(x(:,chunkx{i}),varargin{:});
        for k = 1:numel(thres)
            h{k}(chunkx{i},chunkx{i}) = coef >= thres{k}(1) & coef < thres{k}(2);
        end
    end

end

% for k = 1:numel(thres)
%     h{k} = sparse(h{k});
% end


function maxthr = shiChkMem(px,py,nthr)
if isempty(py)
    py = px;
end
true(px,py,1);
maxthr = 1;
for i = 1:nthr
    try
        true(px,py,i+1);
        maxthr = i;
        clear test;
    catch
        break;
    end
end
