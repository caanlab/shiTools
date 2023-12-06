function [SpmFile,matlabbatch] = shiSpmEstSpm(Dir,Img,Multicond,RegressorTxt,OtherRegressor,OtherRegressorName,TR,MicroTimeResol,MicroTimeOnset,HighPassFilter,inclTimeDeriv,inclDispDeriv,GlobalScaling,Mask,SerialCorr)

% performs first-level estimation
%
% Dir                   - result directory
% Img                   - input images
% Multicond             - multiple condition file name
% RegressorTxt          - covariate file name(s) (e.g. motion param)
% OtherRegressor        - additional covariates (numeric matrix)
% OtherRegressorName    - additional covariate names
% TR                    - in sec (e.g. 2)
% MicroTimeResol        - microtime resolution (e.g. 32 or 2000)
% MicroTimeOnset        - microtime onset (e.g. 16 or 990)
% HighPassFilter        - in sec (e.g. 128)
% inclTimeDeriv         - include time derivatives? (def: false)
% inclDispDeriv         - include dispersion derivatives? (def: false)
% GlobalScaling         - do global scaling? (def: false)
% Mask                  - inclusive mask image file name
% SerialCorr            - 'AR(1)' or 'FAST' (def: 'FAST')

if ~exist('RegressorName','var') || isempty(OtherRegressorName)
    OtherRegressorName = cell(size(OtherRegressor,2),1);
    for i = 1:size(OtherRegressor,2)
        OtherRegressorName{i,1} = sprintf('Cov%02d',i);
    end
end
if ~exist('Mask','var') || isempty(Mask)
    Mask = '';
end
if ~exist('SerialCorr','var') || isempty(SerialCorr)
    SerialCorr = 'FAST';
end

Dir = char(Dir);
Img = cellstr(char(Img));
Multicond = char(Multicond);
RegressorTxt = cellstr(char(RegressorTxt));
OtherRegressorName = cellstr(char(Img));
Mask = char(Mask);



%-----------------------------------------------------------------------
% Job saved on 16-Dec-2019 14:58:42 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.fmri_spec.dir = {shiMkdir(Dir)};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = MicroTimeResol;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = MicroTimeOnset;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = Img;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {Multicond};
for i = 1:size(OtherRegressor,2)
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(i).name = OtherRegressorName{i};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(i).val = OtherRegressor(:,i);
end
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = RegressorTxt;
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = HighPassFilter;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [inclTimeDeriv*1 inclDispDeriv*1];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = shiIf(GlobalScaling,'Scaling','None');
matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
matlabbatch{1}.spm.stats.fmri_spec.mask = {Mask};
matlabbatch{1}.spm.stats.fmri_spec.cvi = SerialCorr;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('serial',matlabbatch);

SpmFile = fullfile(Dir,'SPM.mat');