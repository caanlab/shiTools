function [outImg,TivCsvSaveName,matlabbatch] = shiSpmDartel(c1Img,rc1Img,rc2Img,Seg8Mat,TivCsvSaveName,TemplateName)

% runs DARTEL after segment (dartel create template -> dartel normalize+smooth -> calculate TIV)

c1Img = cellstr(char(c1Img));
rc1Img = cellstr(char(rc1Img));
rc2Img = cellstr(char(rc2Img));
Seg8Mat = cellstr(char(Seg8Mat));

[pth,nme,ext] = shiFileParts(c1Img);
outImg = shiStrConcat(pth,filesep,'smw',nme,ext);

st = shiTime;
if ~exist('TivCsvSaveName','var') || isempty(TivCsvSaveName)
    TivCsvSaveName = fullfile(fileparts(c1Img{1}),['DartelTIV_',st,'.csv']);
end
if ~exist('TemplateName','var') || isempty(TemplateName)
    TemplateName = ['DartelTemplate_',st];
end

mask_ICV = fullfile(fileparts(which('spm.m')),'tpm','mask_ICV.nii');
%-----------------------------------------------------------------------
% Job saved on 19-Feb-2021 17:18:24 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

matlabbatch{1}.spm.tools.dartel.warp.images{1}(1) = rc1Img;
matlabbatch{1}.spm.tools.dartel.warp.images{2}(1) = rc2Img;
matlabbatch{1}.spm.tools.dartel.warp.settings.template = TemplateName;
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;
matlabbatch{2}.spm.tools.dartel.mni_norm.template(1) = cfg_dep('Run Dartel (create Templates): Template (Iteration 6)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{7}));
matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.flowfields(1) = cfg_dep('Run Dartel (create Templates): Flow Fields', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
matlabbatch{2}.spm.tools.dartel.mni_norm.data.subjs.images{1}(1) = c1Img;
matlabbatch{2}.spm.tools.dartel.mni_norm.vox = [1 1 1];
matlabbatch{2}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{2}.spm.tools.dartel.mni_norm.preserve = 1;
matlabbatch{2}.spm.tools.dartel.mni_norm.fwhm = [8 8 8];
matlabbatch{3}.spm.util.tvol.matfiles(1) = Seg8Mat;
matlabbatch{3}.spm.util.tvol.tmax = 3;
matlabbatch{3}.spm.util.tvol.mask = {[mask_ICV,',1']};
matlabbatch{3}.spm.util.tvol.outf = TivCsvSaveName;

spm_jobman('serial',matlabbatch);

