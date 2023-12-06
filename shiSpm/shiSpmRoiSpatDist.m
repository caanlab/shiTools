function DistMat = shiSpmRoiSpatDist(Img1,Img2,Roi,DistType,varargin)

% calculates spatial distance between nifti images at ROIs
%
% shiSpmRoiSpatDist(Img1,Img2,Roi)
% shiSpmRoiSpatDist(Img1,Img2,Roi,Type)
% 
%   Img1     - image(s) of condition 1
%   Img2     - image(s) of condition 2
%   Roi      - one Nifti/MarsBaR ROI file or MNI coordinate+radius (see shiSpmRoiXtr)
%   DistType   - a distance function, or one of the below (or see pdist2)
%      'correlation'   : one minus Pearson correlation, ranging from 0 to 2
%      'spearman'      : one minus Spearman correlation, ranging from 0 to 2
%      'euclidean'     : Euclidean disteance (default)
%      'seuclidean'    : standardized Euclidean distance
%      'absolute'      : absolute distance
%   varargin - fourth and more input to pdist2, if any
%   DistMat  - n_Roi-by-1 cell array of n_Img1-by-n_Img2 matrix
%
% 
%    ###########
% by Zhenhao Shi @ 2018-5-17
%    ###########
% 

if ~exist('DistType','var') || isempty(DistType)
    DistType = 'euclidean';
end
if ~ischar(DistType) && ~isa(DistType,'function_handle')
    error('DistType must be either char or function_handle');
end
if ischar(DistType) && strcmpi(DistType,'absolute')
    DistType = @(x,X)(mean(abs(bsxfun(@minus,x,X)),2));
end

Img1 = cellstr(char(Img1));
Img2 = cellstr(char(Img2));

V1 = spm_vol(char(Img1));
V2 = spm_vol(char(Img2));
[~,XYZ1] = spm_read_vols(V1);
[~,XYZ2] = spm_read_vols(V2);

if ~isequal(XYZ1,XYZ2)
    error('please reslice Image1 and Image2 to match their spaces');
end


X1 = shiSpmRoiXtrVoxelwise(Img1,Roi);
X2 = shiSpmRoiXtrVoxelwise(Img2,Roi);


DistMat = cell(size(X2));

for i = 1:length(X1)
    DistMat{i} = pdist2(X1{i},X2{i},DistType,varargin{:});
end
