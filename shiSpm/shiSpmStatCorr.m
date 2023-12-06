function shiSpmStatCorr(Dir,Img,Vector,Mask,write4D)

% conducts whole-brain correlation analysis and saves R map and Fisher-transformed Z map
%
% shiSpmStatCorr(Dir,Img,Vector)
% 
%   Dir    - string, where results are to be saved
%   Img    - input nifti images
%   Vector - a vector of the same length as Img to correlate voxelwisely
%            with Img
% 
%    ###########
% by Zhenhao Shi @ 2018-1-6
%    ###########
% 

PWD = pwd;

shiMkdir(Dir);
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



matlabbatch{1}.spm.stats.factorial_design.dir = {Dir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = Img;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = Vector;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'Variable';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
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
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [0 1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Neg';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [0 -1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('serial',matlabbatch);

TImg = fullfile(Dir,'spmT_0001.nii');
if ~exist(TImg,'file')
    TImg = fullfile(Dir,'spmT_0001.img');
end
RImg = fullfile(Dir,'corrR.nii');
ZImg = fullfile(Dir,'corrZ.nii');


shiSpmT2R(TImg,RImg,length(Vector));
shiSpmR2Z(RImg,ZImg);


%%
StatType = 'Corr';
Time = shiTime;
if ~exist('write4D','var') || isempty(write4D) || write4D
    shiSpm3dTo4d(Img,fullfile(Dir,'ConImg4d.nii'));
end
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','Mask','StatType','PWD','Time','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','Mask','StatType','PWD','Time');
end
