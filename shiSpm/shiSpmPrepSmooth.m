function outImg = shiSpmPrepSmooth(Img,Fwhm,Prefix)

% spatially smoothes nifti image files with FWHM of 8mm
%
% outImg = shiSpmPrepSmth(Img)
% outImg = shiSpmPrepSmth(Img,Fwhm)
% outImg = shiSpmPrepSmth(Img,Fwhm,Prefix)
%
%   Img         - string or cell array of strings for raw image file names
%   Fwhm        - scalar or 1-by-3 vector of FWHM in mm (default = 8)
%   Prefix      - a char array for output image prefix (default = 's')
%   outImg      - smoothed image file names
% 
%    ###########
% by Zhenhao Shi @ 2017-5-9
%    ###########
% 

warning('obsolete. use shiSpmPreprocSmooth');

Img = cellstr(char(Img));

if ~exist('Fwhm','var') || isempty(Fwhm)
    Fwhm = 8;
end
if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 's';
end

[pth,nme,ext] = shiFileParts(Img);

outImg = cell(size(Img));
for i = 1:length(Img)
    outImg{i} = fullfile(pth{i},[Prefix,nme{i},ext{i}]);
end

for i = 1:length(Img)
    spm_smooth(Img{i},outImg{i},Fwhm)
end


% %-----------------------------------------------------------------------
% % Job configuration created by cfg_util (rev $Rev: 4252 $)
% %-----------------------------------------------------------------------
% matlabbatch{1}.spm.spatial.smooth.data = Img;
% matlabbatch{1}.spm.spatial.smooth.fwhm = Fwhm;
% matlabbatch{1}.spm.spatial.smooth.dtype = 0;
% matlabbatch{1}.spm.spatial.smooth.im = 0;
% matlabbatch{1}.spm.spatial.smooth.prefix = Prefix;
% 
% spm_jobman('serial',matlabbatch);
% 
