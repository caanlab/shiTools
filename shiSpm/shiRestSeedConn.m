function shiRestSeedConn(PathData,PathOut,Roi,RoiName)

% conduts seed-based connectivity analysis on resting-state data
% 
% shiRestSeedConn(Dir,Img,Roi)
% shiRestSeedConn(Dir,Img,Roi,Radius)
% 
%   Dir         - string, path to save the results
%   Img         - cell array of strings, .img/.nii data file names
%   Roi         - string of ROI file name (.img/.nii/), or MNI coordinates,
%                 or coordinates+radius (see shiSpmRoiXtr)
%   RoiName     - string for file names (RoiName_beta.img, RoiName_Z.img)
% 
% 
%    ###########
% by Zhenhao Shi @ 2017-5-13
%    ###########
%


Img = shiFullFileName(fullfile(PathData,'*.img'));
if isempty(Img)
    Img = shiFullFileName(fullfile(PathData,'*.IMG'));
end
if isempty(Img)
    Img = shiFullFileName(fullfile(PathData,'*.nii'));
end
if isempty(Img)
    Img = shiFullFileName(fullfile(PathData,'*.NII'));
end
if isempty(Img) || length(Img)<4
    error('no or too few images found');
end

X=shiSpmRoiXtr(Img,Roi);

tmpDir = shiMkdir(fullfile(PathOut,['temp',num2str(rand())]));

[ZImg,bImg,AnyMiss] = xxx_shiSpmStatCorr(tmpDir,Img,X);

shiSpmImgRename(ZImg,fullfile(PathOut,['Seed_',RoiName,'_Z.img']));
shiSpmImgRename(bImg,fullfile(PathOut,['Seed_',RoiName,'_beta.img']));

rmdir(tmpDir,'s');

if AnyMiss
    Excluded_Image = Img(isnan(X));
    save(fullfile(PathOut,['Seed_',RoiName,'_Excluded_Image.mat']),'Excluded_Image');
end





function [ZImg,bImg,AnyMiss] = xxx_shiSpmStatCorr(tmpDir,Img,Vector)


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


matlabbatch{1}.spm.stats.factorial_design.dir = {tmpDir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = Img;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = Vector;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'Variable';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
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
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('serial',matlabbatch);

TImg = fullfile(tmpDir,'spmT_0001.img');
RImg = fullfile(tmpDir,'corrR.img');
ZImg = fullfile(tmpDir,'corrZ.img');
bImg = fullfile(tmpDir,'con_0001.img');
if ~exist(TImg,'file')
    TImg = fullfile(tmpDir,'spmT_0001.nii');
    RImg = fullfile(tmpDir,'corrR.nii');
    ZImg = fullfile(tmpDir,'corrZ.nii');
    bImg = fullfile(tmpDir,'con_0001.nii');
end
if ~exist(TImg,'file')
    TImg = fullfile(tmpDir,'spmT_0001.IMG');
    RImg = fullfile(tmpDir,'corrR.IMG');
    ZImg = fullfile(tmpDir,'corrZ.IMG');
    bImg = fullfile(tmpDir,'con_0001.IMG');
end
if ~exist(TImg,'file')
    TImg = fullfile(tmpDir,'spmT_0001.NII');
    RImg = fullfile(tmpDir,'corrR.NII');
    ZImg = fullfile(tmpDir,'corrZ.NII');
    bImg = fullfile(tmpDir,'con_0001.NII');
end
if ~exist(bImg,'file') || ~exist(TImg,'file')
    error('cannot find con_0001/spmT_0001 image');
end

shiSpmT2R(TImg,RImg,length(Vector));
shiSpmR2Z(RImg,ZImg);
