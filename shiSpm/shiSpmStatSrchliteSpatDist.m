function shiSpmStatSrchliteSpatDist(Dir,Img,Mask,DistType,Radius,varargin)

% calculates pairwise spatial distance between nifti images using whole-brain searchlights (see pdist)

PWD = pwd;
cd(shiMkdir(Dir));

if ~exist('DistType','var') || isempty(DistType)
    DistType = 'euclidean';
end
if ~ischar(DistType) && ~isa(DistType,'function_handle')
    error('DistType must be either char or function_handle');
end
if ischar(DistType) && strcmpi(DistType,'absolute')
    DistType = @(x,X)(mean(abs(bsxfun(@minus,x,X)),2));
end

Img=cellstr(char(Img));
if ~exist('Mask','var')
    Mask = '';
end
if ~exist('Radius','var')
    Radius = [];
end

spm_vol(char(Img)); % test for error loading images

%%

cnt = 0;
ImgDist_vec = cell(length(Img)*(length(Img)-1)/2,1);
ImgDist_sq = cell(length(Img));
for i = 1:length(Img)-1
    for j = i+1:length(Img)
        cnt = cnt + 1;
        xImg1 = Img{i};
        xImg2 = Img{j};
        xResultImg_name = ['Dist_',num2str(i,'%04d'),'_',num2str(j,'%04d'),'.nii'];
        fprintf('\n%s',xResultImg_name);
        xResultImg = fullfile(Dir,xResultImg_name);
        shiSpmSrchliteSpatDist(xImg1,xImg2,Mask,DistType,Radius,xResultImg,varargin{:});
        ImgDist_vec{cnt} = xResultImg;
        ImgDist_sq{i,j} = xResultImg;
        ImgDist_sq{j,i} = xResultImg;
    end
end
fprintf('\n');

cd(PWD);

%%
StatType = 'SrchliteSpatDist'; %#ok<*NASGU>
if isa(DistType,'function_handle')
    DistType = func2str(DistType);
end
Time = shiTime;
shiSpm3dTo4d(Img,fullfile(Dir,'ConImg4d.nii'));
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','ImgDist_*','StatType','DistType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','ImgDist_*','StatType','DistType','PWD','Time');
end