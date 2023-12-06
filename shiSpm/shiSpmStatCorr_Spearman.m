function shiSpmStatCorr_Spearman(Dir,Img,Mask,Vector)

% conducts whole-brain spearman correlation analysis and saves R map and Fisher-transformed Z map
%
% shiSpmStatCorr_Spearman(Dir,Img,Vector)
%
%   Dir    - string, where results are to be saved
%   Img    - input nifti images
%   Vector - a vector of the same length as Img to correlate voxelwisely
%            with Img
%
%    ###########
% by Zhenhao Shi @ 2018-4-20
%    ###########
%

PWD = pwd;
cd(shiMkdir(Dir));

Img=cellstr(char(Img));

if size(Img,1) ~= size(Vector,1)
    error('unmatched observation number');
end

Img_orig = Img;
Vector_orig = Vector; %#ok<*NASGU>
[Img,Vector,AnyMiss] = shi_deNaN(Img,Vector);

if size(Img,1) ~= size(Vector,1)
    error('unmatched observation number');
end

if ~exist('Mask','var') || isempty(Mask)
    Mask = {''};
else
    Mask = cellstr(char(Mask));
    if numel(Mask)~=1
        error('Mask should be either left empty or specified as one inclusive mask image filename');
    end
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
YR(Mask2) = corr(Y_flat(:,Mask2),Vector(:),'type','spearman');

VR = struct('fname',   'spearmanR.nii',...
            'dim',     V(1).dim(1:3),...
            'dt',      [64 spm_platform('bigend')],...
            'mat',     V(1).mat,...
            'descrip', ['shiSpmStatCorr_Spearman:N=',num2str(N)]);
spm_write_vol(VR,YR);


RImg = fullfile(Dir,'spearmanR.nii');
TImg = fullfile(Dir,'spearmanT.nii');
ZImg = fullfile(Dir,'spearmanZ.nii');

shiSpmR2T(RImg,TImg,length(Vector));
shiSpmR2Z(RImg,ZImg);


%%
StatType = 'Spearman';
Time = shiTime;
shiSpm3dTo4d(Img,fullfile(Dir,'ConImg4d.nii'));
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','Mask','StatType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','Mask','StatType','PWD','Time');
end
