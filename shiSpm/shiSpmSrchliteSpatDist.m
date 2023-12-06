function shiSpmSrchliteSpatDist(Img1,Img2,Mask,DistType,Radius,ResultImg,varargin)

% calculates spatial distance between two nifti images using whole-brain searchlights
%
% shiSpmSrchliteSpatDist(Img1,Img2)
% shiSpmSrchliteSpatDist(Img1,Img2,Mask)
% shiSpmSrchliteSpatDist(Img1,Img2,Mask,DistType)
% shiSpmSrchliteSpatDist(Img1,Img2,Mask,DistType,Radius)
% shiSpmSrchliteSpatDist(Img1,Img2,Mask,DistType,Radius,ResultImg)
% 
%   Img1       - ONE image of condition 1
%   Img2       - ONE image of condition 2
%   Mask       - inclusive mask(s)
%   DistType   - a distance function, or one of the below (or see pdist2)
%      'correlation'   : one minus Pearson correlation, ranging from 0 to 2
%      'spearman'      : one minus Spearman correlation, ranging from 0 to 2
%      'euclidean'     : Euclidean disteance (default)
%      'seuclidean'    : standardized Euclidean distance
%      'absolute'      : absolute distance
%   Radius     - Radius of searchlight spheres (default = 5)
%   ResultImg  - a string of nifti file name (default = 'SpatDist_Img1_Img2.nii')
%   varargin   - fourth and more input to pdist2, if any
%
%    ###########
% by Zhenhao Shi @ 2018-5-20
%    ###########
% 

warning('off','stats:pdist2:ConstantPoints');

if ~exist('ResultImg','var') || isempty(ResultImg)
    [~,i1,~] = fileparts(Img1);
    [~,i2,~] = fileparts(Img2);
    ResultImg = ['Dist_',i1,'_',i2,'.nii'];
end
if ~exist('Radius','var') || isempty(Radius)
    Radius = 5;
end
if ~exist('DistType','var') || isempty(DistType)
    DistType = 'euclidean';
end
if ~ischar(DistType) && ~isa(DistType,'function_handle')
    error('DistType must be either char or function_handle');
end
if ischar(DistType) && strcmpi(DistType,'absolute')
    DistType = @(x,X)(mean(abs(bsxfun(@minus,x,X)),2));
end


if ~strcmpi(ResultImg(end-3:end),'.img') && ~strcmpi(ResultImg(end-3:end),'.nii')
    error('must save results in nifti');
end


V1 = spm_vol(Img1);
V2 = spm_vol(Img2);
[Y1,XYZ1] = spm_read_vols(V1);
[Y2,XYZ2] = spm_read_vols(V2);

if ~isequal(XYZ1,XYZ2) || ~isequal(V1.mat,V2.mat)
    error('please reslice Image1 and Image2 to match their spaces');
end

ShiftVec = shi3dSphereCoord(Radius,V1.mat);
omY1 = nan(size(ShiftVec,2),numel(Y1));
omY2 = omY1;
for i = 1:size(ShiftVec,2)
    shiftY1 = shiMatShift(Y1,NaN,ShiftVec(:,i));
    shiftY2 = shiMatShift(Y2,NaN,ShiftVec(:,i));
    omY1(i,:) = shiftY1(:);
    omY2(i,:) = shiftY2(:);
end

if ~exist('Mask','var') || isempty(Mask)
    Ym = true(size(Y1));
else
    [~,Ym] = shiSpmMaskRead(Img1,Mask,'incl');
end
if any(isnan(Y1(:)))
    mskY1 = ~isnan(Y1);
else
    mskY1 = Y1~=0;
end
if any(isnan(Y2(:)))
    mskY2 = ~isnan(Y2);
else
    mskY2 = Y2~=0;
end
Ym = Ym & mskY1 & mskY2;
Ym(isnan(sum(omY1))) = false;
Ym(isnan(sum(omY2))) = false;


Y3 = nan(size(Y1));
xxx = find(Ym(:));


for i = 1:length(xxx)

    xY1 = omY1(:,xxx(i));
    xY2 = omY2(:,xxx(i));

%     try
        Y3(xxx(i)) = pdist2(xY1(:)',xY2(:)',DistType,varargin{:});
%     catch
%         fprintf(2,' [%g %g %g]: error, NaN returned\n       ',XYZ1(:,xxx(i)));
%         Y3(xxx(i)) = NaN;
%     end
end



if isa(DistType,'function_handle')
    DistType = func2str(DistType);
end

V3 = struct('fname',   ResultImg,...
            'dim',     V1.dim,...
            'dt',      [64 spm_platform('bigend')],...
            'mat',     V1.mat,...
            'descrip', ['shiSpmSrchliteSpatDist(',DistType,')']);

V3.fname = ResultImg;
V3.private.dat.fname = ResultImg;
spm_write_vol(V3,Y3);


warning('on','stats:pdist2:ConstantPoints');

fprintf('\n');





function Coord = shi3dSphereCoord(Radius,AffineMatrix)
A = AffineMatrix(1:3,1:3);
d = 1;
NoMore = false;
Coord = [0;0;0];
while ~NoMore
    d = d+2;
    xCoord = shi3dBoxSurfaceCoord(d);
    xCoord = xCoord(:,sum((A*xCoord).^2)<=Radius^2);
    if isempty(xCoord)
        NoMore = true;
    else
        Coord = [Coord,xCoord];
    end
end


function Coord = shi3dBoxSurfaceCoord(d)
if mod(d,2)~=1 || d<=0
    error('d must be positive odd number');
end
if d == 1
    Coord = [0;0;0];
    return;
end
[iAll,jAll,kAll]=ind2sub([d,d,d],1:d^3);
CoordAll=[iAll;jAll;kAll]-(d+1)/2;
[iInner,jInner,kInner]=ind2sub([d-2,d-2,d-2],1:(d-2)^3);
CoordInner=[iInner;jInner;kInner]-(d-1)/2;
Coord = setdiff(CoordAll',CoordInner','rows')';