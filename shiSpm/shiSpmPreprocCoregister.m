function [outImg,matlabbatch] = shiSpmPreprocCoregister(Img,ImgTemplate,ImgOther,doWrite,Prefix,existAction)

% performs preprocessing: coregistration
%
% outImg = shiSpmPreprocCoregister(Img,ImgTemplate)
% outImg = shiSpmPreprocCoregister(Img,ImgTemplate,ImgOther)
% outImg = shiSpmPreprocCoregister(Img,ImgTemplate,ImgOther,doWrite)
% outImg = shiSpmPreprocCoregister(Img,ImgTemplate,ImgOther,doWrite,Prefix)
%
%   Img           - Raw image (usually structural)
%   ImgTemplate   - Reference image (usually mean realigned functional)
%   ImgOther      - other images to be coregistered to template
%   doWrite       - false (default): est only. true: est & write
%   Prefix        - add prefix only when doWrite==true (ignored otherwise)
%   outImg        - coregistered image (usually structural)
%
% Zhenhao Shi, 2019-10-30
%


Img = cellstr(char(Img));
ImgTemplate = cellstr(char(ImgTemplate));
if ~exist('ImgOther','var') || isempty(ImgOther)
    ImgOther = {''};
else
    ImgOther = cellstr(char(ImgOther));
end

if ~exist('doWrite','var') || isempty(doWrite)
    doWrite = false;
end

if doWrite
    if ~exist('Prefix','var') || isempty(Prefix)
        Prefix = 'r';
    end
    [pth,nme,ext] = shiFileParts(Img);
    outImg = shiStrConcat(pth,filesep,Prefix,nme,ext);
    if ~exist('existAction','var') || isempty(existAction)
        existAction = 'ask';
    elseif ~ismember(lower(existAction),{'ask','overwrite'})
        error('input not recognized');
    end
    if exist(outImg{1},'file') && strcmpi(existAction,'ask')
        ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',outImg{1}), '\', '\\'),'s');
        switch lower(ACTION)
            case 'n'
                error('aborted');
            case 'y'
                warning('overwritting...');
            otherwise
                error('input not recognized');
        end
    end
else
    outImg = Img;
end


%-----------------------------------------------------------------------
% Job saved on 30-Oct-2019 14:06:31 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
if ~doWrite
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = ImgTemplate;
    matlabbatch{1}.spm.spatial.coreg.estimate.source = Img;
    matlabbatch{1}.spm.spatial.coreg.estimate.other = ImgOther;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
else
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = ImgTemplate;
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = Img;
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = ImgOther;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = Prefix;
end

spm_jobman('serial',matlabbatch);
