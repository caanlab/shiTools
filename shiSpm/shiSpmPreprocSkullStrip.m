function [outImg,mskImg] = shiSpmPreprocSkullStrip(Img,c1Img,c2Img,c3Img,Expr,FillHole,Prefix,existAction)

% performs skull stripping by applying tissue masks to functional images
%
% [outImg,mskImg] = shiSpmPreprocSkullStrip(Img,c1Img,c2Img,c3Img)
% [outImg,mskImg] = shiSpmPreprocSkullStrip(Img,c1Img,c2Img,c3Img,Expr)
% [outImg,mskImg] = shiSpmPreprocSkullStrip(Img,c1Img,c2Img,c3Img,Expr,FillHole)
% [outImg,mskImg] = shiSpmPreprocSkullStrip(Img,c1Img,c2Img,c3Img,Expr,FillHole,Prefix)
%
% Img - raw functional images
% c1Img, c2Img, c3Img - c1, c2 and c3 images from segmentation
% Expr - expression applied to c1, c2 and c3 images
%         e.g. 'i1+i2>0.2', for c1+c2>0.2 (SPM)
%              'i1+i2+i3>0.5', for c1+c2+c3>0.5 (FSL) (default)
% FillHole - filling holes in the mask (def = true)
% Prefix - default: 'b'
% outImg - skull stripped functional images
% mskImg - mask used for skull stripping (MaskDebone_*.nii)
%
% Zhenhao Shi 2019/12/10

[pth,nme,ext] = fileparts(Img{1});
xout = fullfile(pth,[Prefix,nme,ext]);

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
c1Img = cellstr(char(c1Img));
c2Img = cellstr(char(c2Img));
c3Img = cellstr(char(c3Img));
mskImg = fullfile(pth,['MaskDebone_',nme,ext]);

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'b';
end

if ~exist('FillHole','var') || isempty(FillHole)
    FillHole = True;
end

if ~exist('Expr','var') || isempty(Expr)
    Expr = 'i1+i2+i3>0.5';
end

shiSpmImgCalc0([c1Img;c2Img;c3Img],Expr,mskImg);

if FillHole
    V = spm_vol(mskImg);
    Y = spm_read_vols(V);
    Y = imfill(Y,'holes');
    spm_write_vol(V,Y);
end

outImg = shiSpmMask(Img,mskImg,'',0,Prefix);
