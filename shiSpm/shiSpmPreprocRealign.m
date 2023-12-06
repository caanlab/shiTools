function [outImg,outMeanImg,outRpTxt,matlabbatch] = shiSpmPreprocRealign(Img,Prefix,existAction)

% performs preprocessing: realignment
%
% [outImg,outMeanImg,outRp] = shiSpmPreprocRealign(Img)
% [outImg,outMeanImg,outRp] = shiSpmPreprocRealign(Img,Prefix)
%
%   Img             - cellstr, functional images
%   outImg          - realigned image filenames (cellstr)
%   outMeanImg      - mean realigned image filename (str)
%   outRpTxt        - motion parameter filename (str)
%
% Zhenhao Shi, 2019-10-30
%


Img = cellstr(char(Img));
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
if exist(outImg{end},'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',outImg{end}), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end

outMeanImg_old = char(shiStrConcat(pth{1},filesep,'mean',nme{1},ext{1}));
outMeanImg = char(shiStrConcat(pth{1},filesep,'Mean_',nme{1},ext{1}));
outRpTxt_old = char(shiStrConcat(pth{1},filesep,'rp_',nme{1},'.txt'));
outRpTxt_old2 = char(shiStrConcat(pth{1},filesep,'rp2_',nme{1},'.txt'));
outRpTxt = char(shiStrConcat(pth{1},filesep,'Rp_',nme{1},'.txt'));

%-----------------------------------------------------------------------
% Job saved on 30-Oct-2019 14:04:10 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.realign.estwrite.data = {Img}';
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = Prefix;

spm_jobman('serial',matlabbatch);
movefile(outMeanImg_old,outMeanImg);
movefile(outRpTxt_old,outRpTxt_old2);
movefile(outRpTxt_old2,outRpTxt);
