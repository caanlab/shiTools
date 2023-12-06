function ImgNorm = shiSpmGlobalNorm(Img,Mask,GlobalNorm)

% performs global normalisation (needs revision)
%
%   GlobalNorm - one of the below:
%                   'none'     - scaled by grand mean (default, see SPM manual)
%                   'scaling'  - scaled by volumn mean (see SPM manual)
% ImgNorm = norm_Img
% (masking procedure will be revised)


warning('this function needs revision')

if nargin < 3 || isempty(GlobalNorm)
    GlobalNorm = 'none';
end

if nargin < 2 || isempty(Mask)
    useSpmFunc = 1;
else
    useSpmFunc = 0;
end

Img = cellstr(char(Img));


G = nan(length(Img),1);

if useSpmFunc
    for i = 1:length(Img)
        G(i,1) = spm_global(spm_vol(Img{i}));
    end;
    ImgMasked = Img;
else
    [ImgMasked,YMask] = shiSpmMask(Img,Mask,'inclusive');
    for i = 1:length(Img)
        validVox = shiNiftiRead_old(ImgMasked{i}).*(YMask>0);
        G(i,1) = nanmean(validVox(:));
    end;
end


if strcmpi(GlobalNorm,'none')
    G(:) = mean(G);
end

ImgNorm = cell(size(Img));
for i = 1:length(Img)
    [xxpath,xxname,xxext] = fileparts(Img{i});
    ImgNorm{i} = fullfile(xxpath,['norm_',xxname,xxext]);
    Func = ['i1./',num2str(G(i),'%.8g'),'.*100'];
    shiSpmImgCalc0(ImgMasked(i),Func,ImgNorm{i});
end

if ~useSpmFunc
    try
        for i = 1:length(Img)
            delete(ImgMasked{i});
            delete([ImgMasked{i}(1:end-4),'.hdr']);
        end;
    catch
        warning('masked but un-normalized images not fully deleted');
    end
end



