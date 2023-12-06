function X = shiSpmRoiXtr_GlobalNorm(Image,Roi,Radius,GlobalNorm)
 
% returns the mean value of all voxels within an ROI and performs global normalization
% 
% X = shiSpmRoiXtr_GlobalNorm(Image,Roi)
% X = shiSpmRoiXtr_GlobalNorm(Image,Roi,Radius)
% X = shiSpmRoiXtr_GlobalNorm(Image,Roi,Radius,GlobalNorm)
% 
%   Image      - a string or a cell array of strings for nifti image file
%                name(s)
%   Roi        - see shiSpmRoiXtr
%   Radius     - a positive number in millimeter, to specify the radius of
%                spheric ROIs, but only when variable Roi contains 
%                coordinates (default = 5)
%   GlobalNorm - one of the below:
%                   'raw'      - do not perform global normalization
%                   'none'     - scaled by grand mean (default, see SPM manual)
%                   'scaling'  - scaled by volumn mean (see SPM manual)
%   X          - each column corresponds to an ROI and each row to an image
% 
%    ###########
% by Zhenhao Shi @ 2018-1-3
%    ###########
% 

if nargin < 3
    Radius = 5;
end
if nargin < 4
    GlobalNorm = 'none';
end

if ~ischar(GlobalNorm) || ~shiStrIncl(GlobalNorm,{'raw','none','scaling'},false)
    error('specify GlobalNorm as ''raw'', ''none'', or ''scaling''');
end

Image=cellstr(char(Image));
X = shiSpmRoiXtr_old(Image,Roi,Radius);

if strcmpi(GlobalNorm,'raw')
    return;
end

G = nan(numel(Image),1);
fprintf('computing global... 0000');
for i = 1:numel(Image)
    fprintf('\b\b\b\b%.4d',i);
    G(i) = spm_global(spm_vol(Image{i}));
end
fprintf('\n');

if strcmpi(GlobalNorm,'none')
    X = X./mean(G).*100;
elseif strcmpi(GlobalNorm,'scaling')
    X = X./G.*100;
end
    