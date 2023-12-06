function [outImg,matlabbatch] = shiSpmPreprocSmooth(Img,Fwhm,Prefix,existAction)

% performs preprocessing: spatially smoothing
%
% outImg = shiSpmPreprocSmooth(Img)
% outImg = shiSpmPreprocSmooth(Img,Fwhm)
% outImg = shiSpmPreprocSmooth(Img,Fwhm,Prefix)
%
%   Img         - string or cell array of strings for raw image file names
%   Fwhm        - scalar or 1-by-3 vector of FWHM in mm (default = 8)
%   Prefix      - a char array for output image prefix (default = 's')
%   outImg      - smoothed image file names
% 
% Zhenhao Shi, 2019-10-30
% 


Img = cellstr(char(Img));

if ~exist('Fwhm','var') || isempty(Fwhm)
    Fwhm = 8;
end
if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 's';
end

[pth,nme,ext] = shiFileParts(Img);
outImg = shiStrConcat(pth,filesep,Prefix,nme,ext);

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

if isscalar(Fwhm)
    Fwhm = [Fwhm,Fwhm,Fwhm];
end

%-----------------------------------------------------------------------
% Job saved on 30-Oct-2019 14:11:05 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.smooth.data = Img;
matlabbatch{1}.spm.spatial.smooth.fwhm = Fwhm;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = Prefix;

spm_jobman('serial',matlabbatch);
