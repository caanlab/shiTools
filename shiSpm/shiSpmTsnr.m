function outImg = shiSpmTsnr(Img,DetrendOrder,outImg)

% calculates temporal SNR image
%
% Zhenhao Shi 2020-5-13
%

disp('TSNR: reading input...');

Img = cellstr(char(Img));
[pth,nme,ext] = shiFileParts(Img);

if ~exist('outImg','var') || isempty(outImg)
    outImg = fullfile(pth{1},['q',nme{1},ext{1}]);
end
if ~exist('DetrendOrder','var') || isempty(DetrendOrder)
    DetrendOrder = 0;
end

if DetrendOrder > 0
    xxdetrend = @(x)(bsxfun(@plus,mean(x,1),spm_detrend(x,DetrendOrder)));
else
    xxdetrend = @(x)(x);
end


V = spm_vol(char(Img));
Y = spm_read_vols(V);

Y = reshape(Y,prod(V(1).dim(1:3)),length(Img))';

disp('TSNR: calculating tSNR...');

Y = xxdetrend(Y);
Y = mean(Y)./std(Y);
Yo = reshape(Y,V(1).dim);

disp('TSNR: writing...');

    Vo = V(1);
    Vo.fname = char(outImg);
    Vo.descrip = [Vo.descrip, sprintf('; tSNR, DetrendOrder=%g', DetrendOrder)];
    spm_write_vol(Vo,Yo);

disp('TSNR: finished.');
