function shiSpmStatDist(Dir,Img,Vector,DistType,Mask,varargin)

% conducts whole-brain distance analysis for the distance between images and a vector, and returns distance map
%
% shiSpmStatDist(Dir,Img,Vector)
%
%   Dir    - string, where results are to be saved
%   Img    - a cell array input nifti images (or, for the purpose of RSA, cells as distance matrix)
%   Vector - a vector of the same length as Img to compute distance voxelwisely against Img
%            with Img (or, for the purpose of RSA, distance matrix)
%   DistType   - a distance function, or one of the below
%      'correlation'   : one minus Pearson correlation, ranging from 0 to 2
%      'spearman'      : one minus Spearman correlation, ranging from 0 to 2
%      'euclidean'     : Euclidean disteance (default)
%      'seuclidean'    : standardized Euclidean distance
%      'mahalanobis'   : Mahalanobis distance
%      'absolute'      : absolute distance
%   varargin - fourth and more input to pdist, if any
%
%    ###########
% by Zhenhao Shi @ 2018-5-18
%    ###########
%

PWD = pwd;
cd(shiMkdir(Dir));

if iscell(Img) && isequal(Img,Img') % Img is arranged as a symmetric distance matrix
    warning('transforming Img from matrix to array');
    Img_square = Img;
    Img = cell(size(Img,1)*(size(Img,1)-1)/2,1);
    cnt = 0;
    for i = 1:size(Img_square,1)-1
        for j = i+1:size(Img_square,1)
            cnt = cnt + 1;
            Img{cnt} = Img_square{i,j};
        end
    end
end
if isequal(Vector,Vector') % Vector is arranged as a symmetric distance matrix
    warning('transforming Vector from matrix to array');
    Vector = squareform(Vector)';
end

Img=cellstr(char(Img));

if size(Img,1) ~= size(Vector,1)
    error('unmatched observation number');
end

if ~exist('DistType','var') || isempty(DistType)
    DistType = 'euclidean';
end
if ~ischar(DistType) && isa(DistType,'function_handle')
    error('DistType must be either char or function_handle');
end
if ischar(DistType) && strcmpi(DistType,'absolute')
    DistType = @(x,X)(mean(abs(bsxfun(@minus,x,X)),2));
end

if ~exist('Mask','var') || isempty(Mask)
    Mask = {''};
else
    Mask = cellstr(char(Mask));
    if numel(Mask)~=1
        error('Mask should be either left empty or specified as one inclusive mask image filename');
    end
end


Img_orig = Img;
Vector_orig = Vector; %#ok<*NASGU>
[Img,Vector,AnyMiss] = shi_deNaN(Img,Vector);

if size(Img,1) ~= size(Vector,1)
    error('unmatched observation number');
end


if AnyMiss
    warning('\n\n  #############################\n  ##                         ##\n  ##       - WARNING -       ##\n  ##                         ##\n  ##  Missing values found!  ##\n  ##                         ##\n  #############################\n  ##  %s\n\n',Dir);
    save(fullfile(Dir,'_MissingValuePresent.mat'),'Img_orig','Img','Vector_orig','Vector');
end

N = numel(Vector);
P = char(Img);
V = spm_vol(P);
if isequal(Mask,{''})
    Y = spm_read_vols(V);
else
    [~,~,Y] = shiSpmMaskRead(Img,Mask);
end

Y_flat = reshape(Y,numel(Y(:,:,:,1)),size(Y,4))';

IsNan = any(isnan(Y_flat));
IsConst = std(Y_flat)==0;
Mask2 = ~(IsNan | IsConst);

YR = nan(size(Y(:,:,:,1)));
YR(IsConst) = 0;
YR(Mask2) = pdist2(Y_flat(:,Mask2)',Vector(:)',DistType,varargin{:});


%%
if isa(DistType,'function_handle')
    DistType = func2str(DistType);
end
VR = struct('fname',   'Dist.nii',...
            'dim',     V(1).dim(1:3),...
            'dt',      [64 spm_platform('bigend')],...
            'mat',     V(1).mat,...
            'descrip', ['shiSpmStatDist(',DistType,'):N=',num2str(N)]);
spm_write_vol(VR,YR);



%%
StatType = 'Distance';
Time = shiTime;
shiSpm3dTo4d(Img,fullfile(Dir,'ConImg4d.img'));
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','Mask','StatType','DistType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','Mask','StatType','DistType','PWD','Time');
end
