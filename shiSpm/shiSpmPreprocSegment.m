function [outImg_y,outImg_iy,outImg_SaveMore,matlabbatch] = shiSpmPreprocSegment(Img,SaveMore,existAction)

% performs preprocessing: segmentation
%
% [outImg_y,outImg_iy] = shiSpmPreprocSegment(Img)
%
%   Img           - Raw image (usually structural)
%   outImg_y      - deformation field image (native --> MNI)
%   outImg_iy     - deformation field image (native <-- MNI)
%
% Zhenhao Shi, 2019-10-30
%

Img = char(Img);
[pth,nme,ext] = fileparts(Img);
outImg_y = fullfile(pth,['y_',nme,ext]);
outImg_iy = fullfile(pth,['iy_',nme,ext]);

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end
if exist(outImg_y,'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',outImg_y), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end


if ~exist('SaveMore','var') || isempty(SaveMore)
    SaveMore = false;
end

if SaveMore
    SaveMore = 1;
else
    SaveMore = 0;
end

if SaveMore
    outImg_SaveMore = struct(...
        'BiasField',fullfile(pth,['BiasField_',nme,ext]), ...
        'm'    , fullfile(pth,['m'    ,nme,ext]), ...
        'c1'   , fullfile(pth,['c1'   ,nme,ext]), ...
        'rc1'  , fullfile(pth,['rc1'  ,nme,ext]), ...
        'wc1'  , fullfile(pth,['wc1'  ,nme,ext]), ...
        'mwc1' , fullfile(pth,['mwc1' ,nme,ext]), ...
        'c2'   , fullfile(pth,['c2'   ,nme,ext]), ...
        'rc2'  , fullfile(pth,['rc2'  ,nme,ext]), ...
        'wc2'  , fullfile(pth,['wc2'  ,nme,ext]), ...
        'mwc2' , fullfile(pth,['mwc2' ,nme,ext]), ...
        'c3'   , fullfile(pth,['c3'   ,nme,ext]), ...
        'rc3'  , fullfile(pth,['rc3'  ,nme,ext]), ...
        'wc3'  , fullfile(pth,['wc3'  ,nme,ext]), ...
        'mwc3' , fullfile(pth,['mwc3' ,nme,ext]) ...
        );
else
    outImg_SaveMore = '';
end


%-----------------------------------------------------------------------
% Job saved on 30-Oct-2019 14:09:00 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.preproc.channel.vols = {Img};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [SaveMore SaveMore];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[which('TPM.nii'),',1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [SaveMore SaveMore];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [SaveMore SaveMore];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[which('TPM.nii'),',2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [SaveMore SaveMore];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [SaveMore SaveMore];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[which('TPM.nii'),',3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [SaveMore SaveMore];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [SaveMore SaveMore];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[which('TPM.nii'),',4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[which('TPM.nii'),',5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[which('TPM.nii'),',6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];

spm_jobman('serial',matlabbatch);
