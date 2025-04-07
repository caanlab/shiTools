function ImgResid2 = shiSpmStatResid(Dir,Img,Vector,Mask)

% saves residual images
% 
% shiSpmStatResid(Dir,Img,Vector)
%   Dir        - string, where results are to be saved
%   Img        - input nifti images
%   Vector     - a matrix of the same number of rows as Img to regress
%                voxelwisely against Img. Each column is a variable
% 
%    ###########
% by Zhenhao Shi @ 2024-8-1
%    ###########
% 


shiMkdir(Dir);
Img=cellstr(char(Img));

if size(Img,1) ~= size(Vector,1)
    error('unmatched observation number');
end

Img_orig = Img;
Vector_orig = Vector; %#ok<*NASGU>
[Img,Vector] = shi_deNaN(Img,Vector);

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


matlabbatch{1}.spm.stats.factorial_design.dir = {Dir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = Img;


for j = 1:size(Vector,2)
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).c = Vector(:,j);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).cname = sprintf('var%02d',j);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).iCC = 1;
end
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).c = ones(length(Img),1);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).cname = 'Constant';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).iCC = 5;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
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
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

ImgResid = fullfile(Dir,shiStrConcat('Res_',num2str((1:length(Img))','%04d'),'.nii'));

spm_jobman('serial',matlabbatch);

[~,nme,ext] = shiFileParts(Img);

ImgMean = fullfile(Dir,'Resid_mean.nii');
shiSpmImgCalc0(Img,'mean(X)',ImgMean);

ImgResid2 = cell(size(Img));
if numel(unique(nme)) == numel(nme)
    for i = 1:length(Img)
        ImgResid2{i} = fullfile(Dir,[nme{i},'_resid',ext{i}]);
    end
else
    for i = 1:length(Img)
        ImgResid2{i} = fullfile(Dir,['ImgResid_',num2str(i,'%04d'),ext{i}]);
    end
end

for i = 1:length(Img)
    shiSpmImgCalc0({ImgResid{i};ImgMean},'i1+i2',ImgResid2{i});
end

cellfun(@delete,ImgResid);
delete(fullfile(Dir,'mask.nii'));
delete(fullfile(Dir,'ResMS.nii'));
delete(fullfile(Dir,'RPV.nii'));
delete(fullfile(Dir,'beta_0*.nii'));
delete(ImgMean);


