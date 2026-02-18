function [outImg,mskImg] = shiSpmPreprocDebone(Img,Tissue,Expr,FillHole,Prefix,existAction)

% performs skull stripping by applying a tissue-based mask to functional images
%
% [outImg,mskImg] = shiSpmPreprocDebone(Img, {c1Img;c2Img;c3Img}, ...)
% [outImg,mskImg] = shiSpmPreprocDebone(Img, c0Img, ...)
% [outImg,mskImg] = shiSpmPreprocSkullStrip(Img, Mask, ...)
%
% Img - raw functional images
% c0Img, c1Img, c2Img, c3Img - PVE, c1, c2 and c3 images from segmentation
% Expr - expression applied to c0/mask image (e.g., 'i1>0' (default))
%        or, expression for c1, c2 and c3 images
%        e.g. 'i1+i2>0.2', for c1+c2>0.2 (SPM)
%             'i1+i2+i3>0.5', for c1+c2+c3>0.5 (FSL) (default)
% FillHole - filling holes in the mask (def = true)
% Prefix - default: 'b'
% outImg - skull stripped functional images
% mskImg - mask used for skull stripping (MaskDebone_*.nii)
%
% Zhenhao Shi 2025-5-16

[pth,nme,ext] = fileparts(Img{1});
xout = fullfile(pth,[Prefix,nme,ext]);

Tissue = cellstr(char(Tissue));
if length(Tissue) == 3
    Method = 'c123';
elseif isscalar(Tissue)
    Method = 'c0';
else
    error('second input must be either c1,c2,c3 images or one c0/mask image');
end

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end
if exist(xout,'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',xout), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end


Img = cellstr(char(Img));
[pth,nme,ext] = fileparts(Img{1});
mskImg = fullfile(pth,['MaskDebone_',nme,ext]);

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'b';
end

if ~exist('FillHole','var') || isempty(FillHole)
    FillHole = true;
end

if ~exist('Expr','var') || isempty(Expr)
    switch Method
        case 'c123'
            Expr = 'i1+i2+i3>0.5';
        case 'c0'
            Expr = 'i1>0';
    end
end

shiSpmImgCalc0(Tissue,Expr,mskImg);

if FillHole
    V = spm_vol(mskImg);
    Y = spm_read_vols(V);
    Y = imfill(Y,'holes');
    spm_write_vol(V,Y);
end

outImg = shiSpmMask(Img,mskImg,'',0,Prefix);
