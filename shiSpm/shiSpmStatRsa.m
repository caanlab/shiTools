function shiSpmStatRsa(Dir,Img,BehavDistMat,Mask,DistType_ImgImg,DistType_ImgBehav,Radius,VarArgIn_ImgImg,VarArgIn_ImgBehav)

% runs representational similarity analysis in a pipeline
%

PWD = pwd;

if ~exist('VarArgIn_ImgImg','var') || isempty(VarArgIn_ImgImg)
    VarArgIn_ImgImg = {};
end
if ~exist('VarArgIn_ImgBehav','var') || isempty(VarArgIn_ImgBehav)
    VarArgIn_ImgBehav = {};
end
if ~exist('Radius','var')
    Radius = '';
end
if ~exist('Mask','var')
    Mask = '';
end
if ~exist('DistType_ImgImg','var')
    DistType_ImgImg = '';
end
if ~exist('DistType_ImgBehav','var')
    DistType_ImgBehav = '';
end



fprintf('\n## RSA: LEVEL 1 ##');

%%%%
shiSpmStatSrchliteSpatDist(shiMkdir(fullfile(Dir,'ImgImgDist')),Img,Mask,DistType_ImgImg,Radius,VarArgIn_ImgImg{:});
%%%%

a = load(fullfile(Dir,'ImgImgDist','StatInfo'),'ImgDist_sq');
cd(fullfile(Dir,'ImgImgDist'));
fprintf('## RSA: LEVEL 2 ##\n\n');

%%%%
shiSpmStatDist(Dir,a.ImgDist_sq,BehavDistMat,DistType_ImgBehav,VarArgIn_ImgBehav{:});
%%%%

cd(PWD);

