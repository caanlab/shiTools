function matlabbatch = shiSpmStatOneSampleT(Dir,Img,Cov,CovName,Mask,write4D)

% conducts whole-brain one-sample t-test
% 
% shiSpmStatOneSampleT(Dir,Img)
% shiSpmStatOneSampleT(Dir,Img,Cov,CovName)
%   Dir     - string, where results are to be saved
%   Img     - input nifti images
%   Cov     - a matrix of covariate(s) with the same number of rows as Img
%             Each column is a variable (default = [])
%   CovName - cell array of strings for covariate names (default = '')
% 
% 
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########
% 

PWD = pwd;

shiMkdir(Dir);
Img = cellstr(char(Img));

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
    if size(Cov,1) ~= size(Img,1) || size(Cov,2) ~= numel(CovName)
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

matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = Img;
%%

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
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Pos';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Neg';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('serial',matlabbatch);



%%
StatType = 'OneSampleT'; %#ok<*NASGU>
Time = shiTime;
if ~exist('write4D','var') || isempty(write4D) || ~write4D
else
    shiSpm3dTo4d(Img,fullfile(Dir,'ConImg4d.nii'));
end
% Img = 'ConImg4d.nii';
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Cov','CovName','Mask','StatType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Cov','CovName','Mask','StatType','PWD','Time');
end

