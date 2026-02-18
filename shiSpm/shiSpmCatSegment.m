function [outImg_y,outImg_iy,out_SaveMore,matlabbatch] = shiSpmCatSegment(Img,SaveMore,OutputLikeSpm,existAction)

% performs CAT12 segmentation (e.g. for VBM)
%
% [outImg_y,outImg_iy] = shiSpmCatSegment(Img)
%
%   Img           - Raw image (usually structural)
%   outImg_y      - deformation field image (native --> MNI)
%   outImg_iy     - deformation field image (native <-- MNI)
%
% Zhenhao Shi, 2025-5-15
%

Img = char(Img);
[pth,nme,ext] = fileparts(Img);
outImg_y = fullfile(pth,'mri',['y_',nme,ext]);
outImg_iy = fullfile(pth,'mri',['iy_',nme,ext]);

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


if ~exist('OutputLikeSpm','var') || isempty(OutputLikeSpm)
    SaveMore = false;
end

if ~exist('SaveMore','var') || isempty(SaveMore)
    SaveMore = false;
end

if SaveMore
    SaveMore = 1;
else
    SaveMore = 0;
end


IQR = fullfile(pth,'report',['IQR_',nme,'.txt']);
TIV = fullfile(pth,'report',['TIV_',nme,'.txt']);


matlabbatch{1}.spm.tools.cat.estwrite.data = {Img};
matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {which('TPM.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasacc = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.optimal = [1 0.3];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.setCOM = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.affmod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASmyostr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.WMHC = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {which('Template_0_GS.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI = struct([]);
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = SaveMore;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;

spm_jobman('serial',matlabbatch);

if SaveMore
    out_SaveMore = struct(...
        'p0'   , fullfile(pth,'mri',['p0'   ,nme,          ext]), ... c0
        'wp0'  , fullfile(pth,'mri',['wp0'  ,nme,          ext]), ... wc0
        'p1'   , fullfile(pth,'mri',['p1'   ,nme,          ext]), ... c1
        'rp1'  , fullfile(pth,'mri',['rp1'  ,nme,'_rigid', ext]), ... rc1
        'mwp1' , fullfile(pth,'mri',['mwp1' ,nme,          ext]), ... mwc1
        'p2'   , fullfile(pth,'mri',['p2'   ,nme,          ext]), ... c2
        'rp2'  , fullfile(pth,'mri',['rp2'  ,nme,'_rigid', ext]), ... rc2
        'mwp2' , fullfile(pth,'mri',['mwp2' ,nme,          ext]), ... mwc2
        'p3'   , fullfile(pth,'mri',['p3'   ,nme,          ext]), ... c3
        'rp3'  , fullfile(pth,'mri',['rp3'  ,nme,'_rigid', ext]), ... rc3
        'mwp3' , fullfile(pth,'mri',['mwp3' ,nme,          ext]), ... mwc3
        'wj'   , fullfile(pth,'mri',['wj_'  ,nme,          ext]), ... wj
        'wm'   , fullfile(pth,'mri',['wm'   ,nme,          ext]), ... wm
        'cat_mat'   , fullfile(pth,'report',['cat_',nme,'.mat']), ...
        'cat_xml'   , fullfile(pth,'report',['cat_',nme,'.xml']), ...
        'catreport' , fullfile(pth,'report',['catreport_',nme,'.pdf']), ...
        'IQR'  , IQR ,...
        'TIV'  , TIV ...
        );
    if ~exist(out_SaveMore.rp1,'file'), out_SaveMore.rp1 = fullfile(pth,'mri',['rp1',nme,ext]); end
    if ~exist(out_SaveMore.rp2,'file'), out_SaveMore.rp2 = fullfile(pth,'mri',['rp2',nme,ext]); end
    if ~exist(out_SaveMore.rp3,'file'), out_SaveMore.rp3 = fullfile(pth,'mri',['rp3',nme,ext]); end
    S = load(out_SaveMore.cat_mat);
    writematrix(S.S.qualityratings.IQR,IQR,'Delimiter','\t');
    writematrix([S.S.subjectmeasures.vol_TIV,S.S.subjectmeasures.vol_abs_CGW([2,3,1,4])],TIV,'Delimiter','\t');
    fldnm = fieldnames(out_SaveMore);
    for f = 1:length(fldnm)
        if ~exist(out_SaveMore.(fldnm{f}),'file')
            out_SaveMore.(fldnm{f}) = '';
        end
    end
else
    out_SaveMore = '';
end


if OutputLikeSpm
    if SaveMore
        out_SaveMore2 = struct(...
            'c0'        , fullfile(pth,['c0'        ,nme,   ext]), ...
            'wc0'       , fullfile(pth,['wc0'       ,nme,   ext]), ...
            'c1'        , fullfile(pth,['c1'        ,nme,   ext]), ...
            'rc1'       , fullfile(pth,['rc1'       ,nme,   ext]), ...
            'mwc1'      , fullfile(pth,['mwc1'      ,nme,   ext]), ...
            'c2'        , fullfile(pth,['c2'        ,nme,   ext]), ...
            'rc2'       , fullfile(pth,['rc2'       ,nme,   ext]), ...
            'mwc2'      , fullfile(pth,['mwc2'      ,nme,   ext]), ...
            'c3'        , fullfile(pth,['c3'        ,nme,   ext]), ...
            'rc3'       , fullfile(pth,['rc3'       ,nme,   ext]), ...
            'mwc3'      , fullfile(pth,['mwc3'      ,nme,   ext]), ...
            'wj'        , fullfile(pth,['wj'        ,nme,   ext]), ...
            'wm'        , fullfile(pth,['wm'        ,nme,   ext]), ...
            'cat_mat'   , fullfile(pth,['CatMat_'   ,nme,'.mat']), ...
            'cat_xml'   , fullfile(pth,['CatXml_'   ,nme,'.xml']), ...
            'catreport' , fullfile(pth,['CatReport_',nme,'.pdf']), ...
            'IQR'       , fullfile(pth,['CatIqr_'   ,nme,'.txt']), ...
            'TIV'       , fullfile(pth,['CatTiv_'   ,nme,'.txt']) ...
            );
        if ~isempty(out_SaveMore.('p0'       )), movefile(out_SaveMore.('p0'       ), out_SaveMore2.('c0'       )); else out_SaveMore2.('c0'       ) = ''; end %#ok<*SEPEX>
        if ~isempty(out_SaveMore.('wp0'      )), movefile(out_SaveMore.('wp0'      ), out_SaveMore2.('wc0'      )); else out_SaveMore2.('wc0'      ) = ''; end
        if ~isempty(out_SaveMore.('p1'       )), movefile(out_SaveMore.('p1'       ), out_SaveMore2.('c1'       )); else out_SaveMore2.('c1'       ) = ''; end
        if ~isempty(out_SaveMore.('rp1'      )), movefile(out_SaveMore.('rp1'      ), out_SaveMore2.('rc1'      )); else out_SaveMore2.('rc1'      ) = ''; end
        if ~isempty(out_SaveMore.('mwp1'     )), movefile(out_SaveMore.('mwp1'     ), out_SaveMore2.('mwc1'     )); else out_SaveMore2.('mwc1'     ) = ''; end
        if ~isempty(out_SaveMore.('p2'       )), movefile(out_SaveMore.('p2'       ), out_SaveMore2.('c2'       )); else out_SaveMore2.('c2'       ) = ''; end
        if ~isempty(out_SaveMore.('rp2'      )), movefile(out_SaveMore.('rp2'      ), out_SaveMore2.('rc2'      )); else out_SaveMore2.('rc2'      ) = ''; end
        if ~isempty(out_SaveMore.('mwp2'     )), movefile(out_SaveMore.('mwp2'     ), out_SaveMore2.('mwc2'     )); else out_SaveMore2.('mwc2'     ) = ''; end
        if ~isempty(out_SaveMore.('p3'       )), movefile(out_SaveMore.('p3'       ), out_SaveMore2.('c3'       )); else out_SaveMore2.('c3'       ) = ''; end
        if ~isempty(out_SaveMore.('rp3'      )), movefile(out_SaveMore.('rp3'      ), out_SaveMore2.('rc3'      )); else out_SaveMore2.('rc3'      ) = ''; end
        if ~isempty(out_SaveMore.('mwp3'     )), movefile(out_SaveMore.('mwp3'     ), out_SaveMore2.('mwc3'     )); else out_SaveMore2.('mwc3'     ) = ''; end
        if ~isempty(out_SaveMore.('wj'       )), movefile(out_SaveMore.('wj'       ), out_SaveMore2.('wj'       )); else out_SaveMore2.('wj'       ) = ''; end
        if ~isempty(out_SaveMore.('wm'       )), movefile(out_SaveMore.('wm'       ), out_SaveMore2.('wm'       )); else out_SaveMore2.('wm'       ) = ''; end
        if ~isempty(out_SaveMore.('cat_mat'  )), movefile(out_SaveMore.('cat_mat'  ), out_SaveMore2.('cat_mat'  )); else out_SaveMore2.('cat_mat'  ) = ''; end
        if ~isempty(out_SaveMore.('cat_xml'  )), movefile(out_SaveMore.('cat_xml'  ), out_SaveMore2.('cat_xml'  )); else out_SaveMore2.('cat_xml'  ) = ''; end
        if ~isempty(out_SaveMore.('catreport')), movefile(out_SaveMore.('catreport'), out_SaveMore2.('catreport')); else out_SaveMore2.('catreport') = ''; end
        if ~isempty(out_SaveMore.('IQR'      )), movefile(out_SaveMore.('IQR'      ), out_SaveMore2.('IQR'      )); else out_SaveMore2.('IQR'      ) = ''; end
        if ~isempty(out_SaveMore.('TIV'      )), movefile(out_SaveMore.('TIV'      ), out_SaveMore2.('TIV'      )); else out_SaveMore2.('TIV'      ) = ''; end
        out_SaveMore = out_SaveMore2;
    end
    outImg_y2 = fullfile(pth,['y_',nme,ext]);
    movefile(outImg_y, outImg_y2);
    outImg_iy2 = fullfile(pth,['iy_',nme,ext]);
    movefile(outImg_iy, outImg_iy2);
    outImg_y = outImg_y2;
    outImg_iy = outImg_iy2;
    rmdir(fullfile(pth,'mri'));
    rmdir(fullfile(pth,'report'),'s');
end
