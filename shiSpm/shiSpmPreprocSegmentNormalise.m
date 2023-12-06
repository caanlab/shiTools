function [outImg,matlabbatch] = shiSpmPreprocSegmentNormalise(ImgStruct,ImgFunc,VoxelSize,BoundingBox,Interpolation,Prefix,existAction)

% performs two preprocessing steps: segmentation, normalization (segmentation outputs will be moved to functional image folder)
%
% Zhenhao Shi, 2022-9-28
%


if ~exist('VoxelSize','var') || isempty(VoxelSize)
    VoxelSize = [3 3 3];
elseif numel(VoxelSize)==1
    VoxelSize = [VoxelSize VoxelSize VoxelSize];
end

if ~exist('BoundingBox','var') || isempty(BoundingBox)
    BoundingBox = [-78 -112 -70;78 76 85];
elseif ischar(BoundingBox) && strcmpi(BoundingBox,'deform')
    BoundingBox = spm_get_bbox(spm_vol(ImgDeform_y));
end

if ~exist('Interpolation','var') || isempty(Interpolation)
    Interpolation = 4;
end

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'w';
end

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end

[ImgY,ImgIY,matlabbatch_1] = shiSpmPreprocSegment(ImgStruct,false,existAction);
[outImg,matlabbatch_2] = shiSpmPreprocNormalise(ImgFunc,ImgY,VoxelSize,BoundingBox,Interpolation,Prefix,existAction);

if ~isequal(fileparts(ImgY),fileparts(outImg{1}))
    try
        movefile(ImgY,fileparts(outImg{1}));
    catch
        warning('unable to move %s',ImgY);
    end
end
if ~isequal(fileparts(ImgIY),fileparts(outImg{1}))
    try
        movefile(ImgIY,fileparts(outImg{1}));
    catch
        warning('unable to move %s',ImgIY);
    end
end
matlabbatch = [matlabbatch_1;matlabbatch_2];
