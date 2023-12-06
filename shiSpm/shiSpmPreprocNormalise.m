function [outImg,matlabbatch] = shiSpmPreprocNormalise(Img,ImgDeform_y,VoxelSize,BoundingBox,Interpolation,Prefix,existAction)

% performs preprocessing: segmentation
%
% outImg = shiSpmPreprocNormalise(Img,ImgDeform_y)
% outImg = shiSpmPreprocNormalise(Img,ImgDeform_y,VoxelSize)
% outImg = shiSpmPreprocNormalise(Img,ImgDeform_y,VoxelSize,BoundingBox)
% outImg = shiSpmPreprocNormalise(Img,ImgDeform_y,VoxelSize,BoundingBox,Interpolation)
% outImg = shiSpmPreprocNormalise(Img,ImgDeform_y,VoxelSize,BoundingBox,Interpolation,Prefix)
% 
%   Img                     - string or cell array of strings for raw image file names
%   ImgDeform_y             - y_*.nii or iy_*.nii
%   VoxelSize               - scalar or 1-by-3 vector of output voxel size in mm (default = [3 3 3])
%   BoundingBox             - [-78 -112 -70;78 76 85] (default); or 'Deform' to use the same as ImgDeform_y; or specify 2*3 matrix
%   Interpolation           - 0: Nearest neighbour
%                             1: Trilinear
%                             2: 2nd Degree B-spline
%                             3: 3rd Degree B-Spline 
%                             4: 4th Degree B-Spline (default)
%                             5: 5th Degree B-Spline 
%                             6: 6th Degree B-Spline 
%                             7: 7th Degree B-Spline 
%
% Zhenhao Shi, 2020-5-18
%

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'w';
end

Img = cellstr(char(Img));
ImgDeform_y = char(ImgDeform_y);

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


if ~exist('VoxelSize','var') || isempty(VoxelSize)
    VoxelSize = [3 3 3];
elseif numel(VoxelSize)==1
    VoxelSize = [VoxelSize VoxelSize VoxelSize];
end

if ~exist('Interpolation','var') || isempty(Interpolation)
    Interpolation = 4;
end
if ~exist('BoundingBox','var') || isempty(BoundingBox)
    BoundingBox = [-78 -112 -70;78 76 85];
elseif ischar(BoundingBox) && strcmpi(BoundingBox,'deform')
    BoundingBox = spm_get_bbox(spm_vol(ImgDeform_y));
end



%-----------------------------------------------------------------------
% Job saved on 30-Oct-2019 14:10:29 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {ImgDeform_y};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = Img;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BoundingBox;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = VoxelSize;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = Interpolation;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = Prefix;

spm_jobman('serial',matlabbatch);