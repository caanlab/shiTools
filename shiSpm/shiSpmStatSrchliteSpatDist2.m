function shiSpmStatSrchliteSpatDist2(Dir,ImgA,ImgB,Mask,DistType,Radius,varargin)

% calculates pairwise spatial distance between two groups of nifti images using whole-brain searchlights (see pdist2)

PWD = pwd;
cd(shiMkdir(Dir));

if ~exist('DistType','var') || isempty(DistType)
    DistType = 'euclidean';
end
if ~ischar(DistType) && isa(DistType,'function_handle')
    error('DistType must be either char or function_handle');
end
if ischar(DistType) && strcmpi(DistType,'absolute')
    DistType = @(x,X)(mean(abs(bsxfun(@minus,x,X)),2));
end

ImgA=cellstr(char(ImgA));
ImgB=cellstr(char(ImgB));
if ~exist('Mask','var')
    Mask = '';
end
if ~exist('Radius','var')
    Radius = [];
end

spm_vol(char([ImgA;ImgB])); % test for error loading images

%%
ImgDist = cell(length(ImgA),length(ImgB));
for i = 1:length(ImgA)
    for j = 1:length(ImgB)
        xImg1 = ImgA{i};
        xImg2 = ImgB{j};
        xResultImg_name = ['Dist_A',num2str(i,'%04d'),'_B',num2str(j,'%04d'),'.nii'];
        xResultImg = fullfile(Dir,xResultImg_name);
        shiSpmSrchliteSpatDist(xImg1,xImg2,Mask,DistType,Radius,xResultImg,varargin{:});
        ImgDist{i,j} = xResultImg;
    end;
end;

cd(PWD);

%%
StatType = 'SrchliteSpatDist'; %#ok<*NASGU>
if isa(DistType,'function_handle')
    DistType = func2str(DistType);
end
Time = shiTime;
shiSpm3dTo4d(ImgA,fullfile(Dir,'ConImg4d_A.nii'));
shiSpm3dTo4d(ImgB,fullfile(Dir,'ConImg4d_B.nii'));
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','ImgA','ImgB','ImgDist','StatType','DistType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','ImgA','ImgB','ImgDist','StatType','DistType','PWD','Time');
end

