function [outImg,outSpmMat,matlabbatch] = shiSpmPreprocRegressOut(Img,CovData,CovTxt,CovRoi,Prefix,existAction)

% performs band-pass filtering
%
% outImg = shiSpmPreprocRegressOut(Img)
% outImg = shiSpmPreprocRegressOut(Img,CovData)
% outImg = shiSpmPreprocRegressOut(Img,CovData,CovTx)
% outImg = shiSpmPreprocRegressOut(Img,CovData,CovTxt,CovRoi)
% outImg = shiSpmPreprocRegressOut(Img,CovData,CovTxt,CovRoi,Prefix)
% 
%   Img         - raw images
%   CovData     - matrix containing variables to be regressed out (one column per variable) (def: [] (no variable))
%   CovTxt      - text files containing variables to be regressed out (one column per variable) (def: '' (no file))
%   CovRoi      - ROIs whose signals to be regressed out (def: '' (no Roi))
%
% Zhenhao Shi, 2025-3-20
%

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'v';
end
if ~exist('CovRoi','var') || isempty(CovRoi)
    CovRoi = '';
else
    CovRoi = shiSpmRoiFormat(CovRoi);
end
if ~exist('CovTxt','var') || isempty(CovTxt)
    CovTxt = '';
else
    CovTxt = cellstr(char(CovTxt));
end
if ~exist('CovData','var') || isempty(CovData)
    CovData = [];
end


Img = cellstr(char(Img));
[pth,nme,ext] = shiFileParts(Img);
outImg = shiStrConcat(pth,filesep,Prefix,nme,ext);
DirReg = shiMkdir(fullfile(pth{1},[Prefix,nme{1},'_RegressingOut']));

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end
if exist(outImg{1},'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',outImg{1}), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end


% confoundings
Cov_1 = CovData;
Cov_2 = [];
for i = 1:length(CovTxt)
    Cov_2 = [Cov_2,importdata(CovTxt{i})]; %#ok<*AGROW>
end
Cov_3 = shiSpmRoiXtr(Img,CovRoi);
CovAll = [Cov_1, Cov_2, Cov_3];


CovName = [
    shiIf(size(Cov_1,2)==0, [], shiStrConcat('CovDat',1:size(Cov_1,2)));
    shiIf(size(Cov_2,2)==0, [], shiStrConcat('CovTxt',1:size(Cov_2,2)));
    shiIf(size(Cov_3,2)==0, [], shiStrConcat('CovRoi',1:size(Cov_3,2)));
    ];


%-----------------------------------------------------------------------
% Job saved on 31-Oct-2019 11:09:01 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.factorial_design.dir = {DirReg};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = Img;
for i = 1:size(CovAll,2)
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c = CovAll(:,i);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname = CovName{i};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC = 1;
end
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('serial',matlabbatch);

save(fullfile(DirReg,'SPM.mat'),'CovData','CovTxt','CovRoi','-append');

fprintf('moving and renaming residual files... \n');
ImgResid = shiStrConcat(DirReg,filesep,'Res_',cellstr(reshape(sprintf('%04d',1:length(outImg)),4,[])'),'.nii');
SpmMat = fullfile(DirReg,'SPM.mat');
for i = 1:length(outImg)
    movefile(ImgResid{i},outImg{i});
end
outSpmMat = fullfile(pth{1},[Prefix,'Spm_',nme{1},'.mat']);
movefile(SpmMat,outSpmMat);

rmdir(DirReg,'s');