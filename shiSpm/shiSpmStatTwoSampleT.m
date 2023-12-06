function matlabbatch = shiSpmStatTwoSampleT(Dir,Img1,Img2,Cov,CovName,Mask,write4D,UnequalVariance)

% conducts whole-brain paired-sample t-test
%
% shiSpmStatTwoSampleT(Dir,Img1,Img2)
% shiSpmStatTwoSampleT(Dir,Img1,Img2,Cov,CovName)
% 
%   Dir    - string, where results are to be saved
%   Img1   - input nifti images for group 1
%   Img2   - input nifti images for group 2
%   Cov     - a matrix of covariate(s) that correspond to [Img1;Img2]
%             Each column is a variable (default = [])
%   CovName - cell array of strings for covariate names (default = '')
% 
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########
% 

PWD = pwd;

if ~exist('UnequalVariance','var') || isempty(UnequalVariance)
    UnequalVariance = true;
end

shiMkdir(Dir);
Img1=cellstr(char(Img1));
Img2=cellstr(char(Img2));

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
    if size(Cov,1) ~= length(Img1)+length(Img2) || size(Cov,2) ~= length(CovName)
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

matlabbatch{1}.spm.stats.factorial_design.dir = {Dir};
%%
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = Img1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = Img2;
%%
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = UnequalVariance;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;

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
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
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
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '01_G1_G2';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = '02_G2_G1';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = '03_G1_pos';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [1 0];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = '04_G1_neg';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.convec = [-1 0];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = '05_G2_pos';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.convec = [0 1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = '06_G2_neg';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.convec = [0 -1];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = '07_both_pos';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.convec = [0.5 0.5];
matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = '08_both_neg';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.convec = [-0.5 -0.5];
matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

spm_jobman('serial',matlabbatch);




%%
StatType = 'TwoSampleT'; %#ok<*NASGU>
if ~exist('write4D','var') || isempty(write4D) || ~write4D
else
    shiSpm3dTo4d(Img1,fullfile(Dir,'ConImg4d_Grp1.nii'));
    shiSpm3dTo4d(Img2,fullfile(Dir,'ConImg4d_Grp2.nii'));
end
% Img1 = 'ConImg4d_Grp1.nii';
% Img2 = 'ConImg4d_Grp2.nii';
Time = shiTime;
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img1','Img2','Cov','CovName','Mask','StatType','PWD','Time','matlabbatch','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img1','Img2','Cov','CovName','Mask','StatType','PWD','Time','matlabbatch');
end

