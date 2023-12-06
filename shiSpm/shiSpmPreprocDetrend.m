function [outImg,MeanImg] = shiSpmPreprocDetrend(Img,Order,doDemean,Prefix,existAction)

% performs detrending (see spm_detrend.m)
%
% Zhenhao Shi 2020-5-13
%

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'd';
end
if ~exist('doDemean','var') || isempty(doDemean)
    doDemean = false;
end
if ~exist('Order','var') || isempty(Order)
    Order = 2;
end

Img = cellstr(char(Img));
[pth,nme,ext] = shiFileParts(Img);
outImg = shiStrConcat(pth,filesep,Prefix,nme,ext);
MeanImg = fullfile(pth{1},['Mean_',nme{1},ext{1}]);

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end
if exist(outImg{1},'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',outImg{1}), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end


%%

if doDemean
    xxdetrend = @(x)(spm_detrend(x,Order));
else
    xxdetrend = @(x)(bsxfun(@plus,mean(x,1),spm_detrend(x,Order)));
end


V = spm_vol(char(Img));
Y = spm_read_vols(V);

Yo = nan(size(Y));

Y = reshape(Y,prod(V(1).dim(1:3)),length(Img))';

disp('DETRENDING in progress...');

Y = xxdetrend(Y);

for i = 1:length(Img)
    Yo(:,:,:,i) = reshape(Y(i,:),V(1).dim);
    Vo = V(i);
    Vo.fname = outImg{i};
    Vo.descrip = [Vo.descrip, sprintf('; detrended, Order=%g, doDemean=%s', Order, doDemean)];
    Vo.dt(1) = 16;
    spm_write_vol(Vo,single(Yo(:,:,:,i)));
end

shiSpmImgCalc0(Img,'mean(X)',MeanImg);

