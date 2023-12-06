function shiSpmStatOneWayWithinAnova_FlexFact(Dir,ImgConditions,Cov,CovName,Mask,write4D)

% conducts whole-brain one-way repeated-measures anova (flexible factorial)
% 
% shiSpmStatOneWayWithinAnova(Dir,ImgGroups)
% shiSpmStatOneWayWithinAnova(Dir,ImgGroups,Cov,CovName)
%   Dir          - string, where results are to be saved
%   ImgGroups    - cells, each containing same number of images for one condition
%   Cov          - a matrix of covariate(s) with the same number of rows as n_Cond*n_Subj (e.g. [sub1_1; sub1_2; sub2_1; sub2_2; sub3_1; sub3_2; ...])
%                  Each column is a variable (default = [])
%   CovName      - cell array of strings for covariate names (default = '')
% 
%    ###########
% by Zhenhao Shi @ 2021-4-8
%    ###########
% 

PWD = pwd;

Dir = shiMkdir(Dir);
if isempty(Dir)
    error('failed to create directory for analysis');
end

for cond = 1:length(ImgConditions)
    ImgConditions{cond} = cellstr(char(ImgConditions{cond}));
end

for cond = 2:length(ImgConditions)
    if length(ImgConditions{cond-1}) ~= length(ImgConditions{cond})
        error('unequal number of images');
    end
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
    if size(Cov,1) ~= size(ImgConditions{1},1) || size(Cov,2) ~= numel(CovName)
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

n_Sub = length(ImgConditions{1});
n_Cond = length(ImgConditions);

ImgAll = cat(2,ImgConditions{:});



matlabbatch{1}.spm.stats.factorial_design.dir = {Dir};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Time';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
for sub = 1:n_Sub
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub).scans = ImgAll(sub,:)';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub).conds = 1:n_Cond;
end
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 2;
if haveCov
    for k = 1:length(CovName)
        matlabbatch{1}.spm.stats.factorial_design.cov(k).c = Cov(:,k);
        matlabbatch{1}.spm.stats.factorial_design.cov(k).cname = CovName{k};
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCFI = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCC = 1;
    end
end
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = Mask;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

cntContrast = 0;
for i = 1:n_Cond-1
    for j = i+1:n_Cond
        cntContrast = cntContrast+1;
        matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.name = sprintf('T%03d_%02d_vs_%02d',cntContrast,i,j);
        xxxCON = zeros(1,n_Cond);xxxCON(i)=1;xxxCON(j)=-1;
        matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.weights = xxxCON;
        matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.sessrep = 'none';
        cntContrast = cntContrast+1;
        matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.name = sprintf('T%03d_%02d_vs_%02d',cntContrast,j,i);
        xxxCON = zeros(1,n_Cond);xxxCON(i)=-1;xxxCON(j)=1;
        matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.weights = xxxCON;
        matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.sessrep = 'none';
    end
end
for i = 1:n_Cond
    cntContrast = cntContrast+1;
    matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.name = sprintf('T%03d_%02d_vs_others',cntContrast,i);
    xxxCON = -ones(1,n_Cond)/(n_Cond-1);xxxCON(i)=1;
    matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.weights = xxxCON;
    matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.sessrep = 'none';
    cntContrast = cntContrast+1;
    matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.name = sprintf('T%03d_others_vs_%02d',cntContrast,i);
    xxxCON = ones(1,n_Cond)/(n_Cond-1);xxxCON(i)=-1;
    matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.weights = xxxCON;
    matlabbatch{3}.spm.stats.con.consess{cntContrast}.tcon.sessrep = 'none';
end

cntContrast = cntContrast+1;
matlabbatch{3}.spm.stats.con.consess{cntContrast}.fcon.name = 'Main_Effect';
matlabbatch{3}.spm.stats.con.consess{cntContrast}.fcon.weights = eye(n_Cond)-1/n_Cond;
matlabbatch{3}.spm.stats.con.consess{cntContrast}.fcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('serial',matlabbatch);

%%
StatType = 'OneWayWithinAnova_FlexFact'; %#ok<*NASGU>
if ~exist('write4D','var') || isempty(write4D) || ~write4D
else
    for cond = 1:n_Cond
        shiSpm3dTo4d(Img1,fullfile(Dir,sprintf('ConImg4d_Cond%02d.nii',cond)));
    end
end
Time = shiTime;
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','ImgConditions','Cov','CovName','Mask','StatType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','ImgConditions','Cov','CovName','Mask','StatType','PWD','Time');
end