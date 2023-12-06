function shiSpmStatIcc31(Dir,ImgMat,r_0)

% conducts intraclass correlation ICC(3,1)
% 
% shiSpmIcc31(Dir,ImgMat)
% shiSpmIcc31(Dir,ImgMat,r_0)
%   Dir     - string, where results are to be saved
%   ImgMat  - input nifti images, arranged as a n-by-k cell matrix, n=#subjects, k=#repetitions
%   r_0     - intraclass correlation under null hypothesis (default:0)
% 
%    ###########
% by Zhenhao Shi @ 2019-5-30
%    ###########
% 

PWD = pwd;

if ~exist('r0','var') || isempty(r_0)
    r_0 = 0;
end

[n,k] = size(ImgMat);


%% SStotal

Img_SStotal = fullfile(shiMkdir(Dir),'SStotal.nii');
shiSpmImgCalc0(ImgMat(:),sprintf('var(X).*%d',n*k-1),Img_SStotal,true,1);


%% MSR & MSC

Img_mean1 = shiStrConcat(Dir,filesep,'mean1_',1:k,'.nii');
Img_mean2 = shiStrConcat(Dir,filesep,'mean2_',1:n,'.nii');
for kk = 1:k
    shiSpmImgCalc0(ImgMat(:,kk),'mean(X)',Img_mean1{kk},true,1);
end
for nn = 1:n
    shiSpmImgCalc0(ImgMat(nn,:)','mean(X)',Img_mean2{nn},true,1);
end

Img_MSR = fullfile(Dir,'MSR.nii');
shiSpmImgCalc0(Img_mean2,sprintf('var(X).*%d',k),Img_MSR,true,1);

Img_MSC = fullfile(Dir,'MSC.nii');
shiSpmImgCalc0(Img_mean1,sprintf('var(X).*%d',n),Img_MSC,true,1);


%% MSE

Img_MSE = fullfile(Dir,'MSE.nii');
shiSpmImgCalc0({Img_SStotal;Img_MSR;Img_MSC},sprintf('(i1-i2.*%d-i3.*%d)./%d',n-1,k-1,(n-1)*(k-1)),Img_MSE,true,1);


%% r

Img_r = fullfile(Dir,'ICC31_r.nii');
shiSpmImgCalc0({Img_MSR;Img_MSE},sprintf('(i1-i2)./(i1+%d.*i2)',k-1),Img_r,true,1);


%% F

Img_F = fullfile(Dir,'ICC31_F.nii');
[Vo,Yo] = shiSpmImgCalc0({Img_MSR;Img_MSE},sprintf('(i1./i2).*%f',(1-r_0)/(1+(k-1)*r_0)),Img_F,false,1);

df1 = n-1;
df2 = (n-1)*(k-1);
Vo.descrip = sprintf('SPM{F_[%g,%g]} - ICC(3,1), n=%d, k=%d, r_0=%g',df1,df2,n,k,r_0);

spm_write_vol(Vo,Yo);


%%

StatType = 'ICC(3,1)'; %#ok<*NASGU>
Time = shiTime;
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','ImgMat','r_0','StatType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','ImgMat','r_0','StatType','PWD','Time');
end

