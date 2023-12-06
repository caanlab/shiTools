function matlabbatch = shiSpmTfce(SpmMat,ContrastIndex,Mask,n_Perm,MultiThreading)

% runs TFCE estimation of specified contrast(s)
%
% Zhenhao Shi
% 2022-01-31

SpmMat = cellstr(char(SpmMat));

if numel(SpmMat)>1
    for i = 1:numel(SpmMat)
        shiSpmTfce(SpmMat{i},ContrastIndex,n_Perm,Mask);
    end
    return
end

SpmMat = SpmMat{1};

if ~exist('ContrastIndex','var') || isempty(ContrastIndex) || any(ContrastIndex<=0)
    ContrastIndex = Inf;
end

if ~exist('Mask','var') || isempty(Mask)
    Mask = '';
else
    Mask = cellstr(char(Mask));
    if numel(Mask)>1
        error('can only work with ONE mask');
    end
    Mask = char(Mask);
end

if ~exist('n_Perm','var') || isempty(n_Perm) || n_Perm<=0
    n_Perm = 5000;
end

if ~exist('MultiThreading','var') || isempty(MultiThreading)
    MultiThreading = shiIf(ispc,false,true);
end


if ~isempty(Mask)

    if ~exist(Mask,'file')
        error('cannot find mask: %s', Mask);
    end

    Mask_spm = fullfile(fileparts(SpmMat),'mask.nii');
    if ~exist(Mask_spm,'file')
        Mask_spm = fullfile(fileparts(SpmMat),'mask.img');
    end

    if ~exist(Mask_spm,'file')
        warning('cannot find the original mask image with the SPM.mat file.');
    else

        Mask2 = fullfile(fileparts(SpmMat),'tmp_mask.nii');
        copyfile(Mask,Mask2);

        Resliced = shiSpmReslice({Mask_spm;Mask2},0,'tmpTFCE_');
        delete(Resliced{1});
        delete(Mask2);

        Mask2_resliced = fullfile(fileparts(SpmMat),'mask_TFCE.nii');
        movefile(Resliced{2},Mask2_resliced);
        
        Mask2_resliced = {Mask2_resliced};

    end

else
    Mask2_resliced = '';
end

matlabbatch{1}.spm.tools.tfce_estimate.spmmat = {SpmMat};
matlabbatch{1}.spm.tools.tfce_estimate.mask = Mask2_resliced;
matlabbatch{1}.spm.tools.tfce_estimate.conspec.titlestr = '';
matlabbatch{1}.spm.tools.tfce_estimate.conspec.contrasts = ContrastIndex;
matlabbatch{1}.spm.tools.tfce_estimate.conspec.n_perm = n_Perm;
matlabbatch{1}.spm.tools.tfce_estimate.nuisance_method = 2;
matlabbatch{1}.spm.tools.tfce_estimate.tbss = 0;
matlabbatch{1}.spm.tools.tfce_estimate.E_weight = 0.5;
matlabbatch{1}.spm.tools.tfce_estimate.singlethreaded = ~MultiThreading;

spm_jobman('serial',matlabbatch);