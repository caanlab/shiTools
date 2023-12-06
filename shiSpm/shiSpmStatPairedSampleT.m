function shiSpmStatPairedSampleT(Dir,Img1,Img2,Cov,CovName,Mask,write4D)

% conducts whole-brain paired-sample t-test
% 
% shiSpmStatPairedSampleT(Dir,Img1,Img2)
% shiSpmStatPairedSampleT(Dir,Img1,Img2,Cov,CovName)
%   Dir     - string, where results are to be saved
%   Img1    - input nifti images for condition 1
%   Img2    - input nifti images for condition 2 (same length as Img1)
%   Cov     - a matrix of covariate(s) with the same number of rows as n_Cond*n_Subj (e.g. [sub1_1; sub1_2; sub2_1; sub2_2; sub3_1; sub3_2; ...])
%             Each column is a variable (default = [])
%   CovName - cell array of strings for covariate names (default = '')
% 
%    ###########
% by Zhenhao Shi @ 2018-8-9
%    ###########
% 

PWD = pwd;

Dir = shiMkdir(Dir);
if isempty(Dir)
    error('failed to create directory for analysis');
end

Img1 = cellstr(char(Img1));
Img2 = cellstr(char(Img2));

if length(Img1) ~= length(Img2)
    error('unequal number of images');
end

if ~exist('Cov','var') || isempty(Cov)
    haveCov = false;
    Cov = [];
    CovName = '';
elseif ~exist('CovName','var') || isempty(CovName)
    haveCov = true;
    CovName = cell(size(Cov,2),1);
    for j = 1:size(Cov,2)
        CovName{j} = ['Cov',num2str(j)];
    end
else
    haveCov = true;
    CovName = cellstr(char(CovName));
end

if haveCov
    if size(Cov,1) ~= size(Img1,1)*2 || size(Cov,2) ~= numel(CovName)
        error('wrong dimension of Cov/CovName');
    end
end

if ~exist('Mask','var') || isempty(Mask)
    Mask = {''};
else
    Mask = cellstr(char(Mask));
    if numel(Mask)~=1
        error('Mask should be either left empty or specified as one inclusive mask image filename');
    end
end

%-----------------------------------------------------------------------
% Job saved on 09-Aug-2018 15:57:07 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.factorial_design.dir = {Dir};
for im = 1:length(Img1)
    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(im).scans = {Img1{im};Img2{im}};
end
matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;

if haveCov
    for k = 1:length(CovName)
        matlabbatch{1}.spm.stats.factorial_design.cov(k).c = Cov(:,k);
        matlabbatch{1}.spm.stats.factorial_design.cov(k).cname = CovName{k};
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCFI = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCC = 1;
    end
else
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = Mask;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Pos';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Neg';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('serial',matlabbatch);


%%
StatType = 'PairedSampleT'; %#ok<*NASGU>
if ~exist('write4D','var') || isempty(write4D) || write4D
else
    shiSpm3dTo4d(Img1,fullfile(Dir,'ConImg4d_Cond1.nii'));
    shiSpm3dTo4d(Img2,fullfile(Dir,'ConImg4d_Cond2.nii'));
end
Time = shiTime;
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img1','Img2','Cov','CovName','Mask','StatType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img1','Img2','Cov','CovName','Mask','StatType','PWD','Time');
end

