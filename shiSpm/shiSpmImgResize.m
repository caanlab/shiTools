function Img_resize = shiSpmImgResize(Img, VoxelSize, BoundingBox, InterpolationOrder)
%
% changes image voxel size and/or bounding box
%
% VoxelSize: 1x3 vector in mm
% BoundingBox: 2x3 matrix in mm
% InterpolationOrder: default = 1, i.e. linear
%
% Based on John Ashburner's reorient.m
% Adapted by Ged Ridgway
% Adapted by Zhenhao Shi


if ~exist('VoxelSize', 'var') || isempty(VoxelSize) || all(isnan(VoxelSize),'all')
    VoxelSize = [];
end

if ~exist('BoundingBox', 'var') || isempty(BoundingBox) || all(isnan(BoundingBox),'all')
    BoundingBox = [];
end

if ~exist('InterpolationOrder', 'var') || isempty(InterpolationOrder)
    InterpolationOrder = 1; % default = linear
end

if iscell(Img) || ~isvector(Img)
    if ischar(Img)
        Img = cellstr(Img);
    end
    Img_resize = cell(size(Img));
    for i = 1:numel(Img)
        Img_resize{i} = shiSpmImgResize(Img{i}, VoxelSize, BoundingBox, InterpolationOrder);
    end
    return;
end


%%

[pth,nme,ext] = fileparts(Img);
Img_resize = fullfile(pth,[nme,'_resize',ext]);

[BoundingBox_orig,VoxelSize_orig] = spm_get_bbox(Img);

if isempty(VoxelSize)
    VoxelSize = VoxelSize_orig;
end

if isempty(BoundingBox)
    BoundingBox = BoundingBox_orig;
end

V = spm_vol(Img);

BBmin = BoundingBox(1,:);
BBmax = BoundingBox(2,:);

% voxel [1 1 1] of output should map to BB mn
% (the combination of matrices below first maps [1 1 1] to [0 0 0])
mat = spm_matrix([BBmin 0 0 0 VoxelSize(:)']) * spm_matrix([-1 -1 -1]);

% voxel-coords of BB mx gives number of voxels required
% (round up if more than a tenth of a voxel over)
dim = ceil(mat \ [BBmax 1]' - 0.1)';

% output image
VO            = V;
VO.fname      = Img_resize;
VO.dim(1:3)   = dim(1:3);
VO.mat        = mat;
VO = spm_create_vol(VO);
for i = 1:dim(3)
    M = mat\V.mat\spm_matrix([0,0,i]); % inv(spm_matrix([0 0 -i])*(inv(VO.mat)\V.mat)); 
    img = spm_slice_vol(V, M, dim(1:2), InterpolationOrder);
    spm_write_plane(VO, img, i);
end

