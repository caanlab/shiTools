function shiSpmPipe_RestPreproc(xStepInd,anat,func,atlas,CustomCov,CustomSpike,xTr,xTa,xSliceOrder,xRefSlice,xDvars_do,xDespike_Method,xInterpolate_Method,xNormalize_VoxelSize)
%
% resting-state data preprocessing pipeline
%     xStepInd % leave empty for 1:177
%     anat
%     func
%     atlas
%     xTr % TR in seconds
%     xTa % TA in seconds
%     xSliceOrder % slice order or slice time (ms) for slice time correction
%     xRefSlice % reference slice (spatial index or slice time) for slice time correction
%     CustomCov = shiIf(~exist('CustomCov','var') || isempty(CustomCov), {}, CustomCov); % cellstr, leave empty or .txt filename WITHOUT the .txt extension. If to be left empty, must be empty cell {}
%     CustomSpike = shiIf(~exist('CustomSpike','var') || isempty(CustomSpike), {}, CustomSpike); % cellstr, leave empty or .txt filename WITHOUT the .txt extension. If to be left empty, must be empty cell {} if ~exist('CustomCov','var') || isempty(CustomCov)
%     xDvars_do = shiIf(~exist('xDvars_do','var')||isempty(xDvars_do), true, xDvars_do);
%     xDespike_Method = shiIf(~exist('xDespike_Method','var')||isempty(xDespike_Method), '', xDespike_Method); % leave empty for default AFNI '3dDespike', or specify 'ArtRepair'
%     xInterpolate_Method = shiIf(~exist('xInterpolate_Method','var')||isempty(xInterpolate_Method), '', xInterpolate_Method); % leave empty for default 'Lomb' (Lomb-Scargle periodogram), or enter method for interp1.m ('linear', 'cubic', 'spline', etc)
%     xNormalize_VoxelSize = shiIf(~exist('xNormalize_VoxelSize','var')||isempty(xNormalize_VoxelSize), [], xNormalize_VoxelSize); % leave empty for default [3 3 3] (mm)

%% variable naming in this script:
%       xxAbcXyz        : temporary variables that are not tracked
%       xAbcXyz         : step-by-step parameters
%       AbcXyz, abc, a  : paths, filenames, prefixes


if ischar(xStepInd) && isequal(lower(xStepInd),'intro')
    xxDiGraph;
    return;
end

if isempty(xStepInd)
    xStepInd = 1:177;
end

%% ========================================================================

disp('...  INITIALIZING ...');

%% USER INPUT  :  DIRECTORY & FILES

assert(all(cellfun(@(x)exist(x,'file'),cellstr(char(anat)))),'cannot find anat image');
assert(all(cellfun(@(x)exist(x,'file'),cellstr(char(func)))),'cannot find func image(s)');

[pth,anat] = fileparts(char(anat)); %%
[xxpth,func] = shiFileParts(cellstr(char(func))); %%

xxVarArgIn = strsplit('atlas,CustomCov,CustomSpike,xTr,xTa,xSliceOrder,xRefSlice,xDvars_do,xDespike_Method,xInterpolate_Method,xNormalize_VoxelSize', ',')';
xxVarExpr = shiStrConcat(xxVarArgIn,' = [];');
for i = 1:length(xxVarArgIn)
    if ~exist(xxVarArgIn{i},'var'), eval(xxVarExpr{i}); end
end


if isempty(atlas)
    atlas = {
        'AAL'
        'AAL_90'
        'AAL3'
        'AAL3_140'
        'HCPMMP'
        'HOA'
        'Power264'
        'Neuromorphometrics_Brain'
        'Schaefer1000Parcel07Net'
        'Schaefer1000Parcel17Net'
        'Yeo7NetLiberal'
        'Yeo17NetLiberal'
        };
end

assert(~isempty(pth),'must include full path in anat and func');
assert(isequal(pth,xxpth{1}),'anat and func images must be in the same folder');
assert(all(ismember(cellstr(char(atlas)),shiSpmAnatLabel)),'atlas must be included in shiMisc and recognized by shiSpmAnatLabel');

CustomCov = shiIf(isempty(CustomCov), {}, CustomCov); % cellstr, leave empty or .txt filename WITHOUT the .txt extension. If to be left empty, must be empty cell {}
CustomSpike = shiIf(isempty(CustomSpike), {}, CustomSpike); % cellstr, leave empty or .txt filename WITHOUT the .txt extension. If to be left empty, must be empty cell {} if isempty(CustomCov)

%% USER INPUT  :  STEP-BY-STEP PREFIXES AND PARAMETERS

% [A] Motion estimation

% [B] Expand motion parameters and get FD-spikes
xMotion_FdOption = []; % leave empty for default 'Fd_Jenkinson';
xMotion_FdSpikeThres = []; % leave empty for default 0.5;
Mot24 = 'Mot24';
Fd = 'Fd';
SpikeFd = 'SpikeFd';

% [C] Adjust for slice timing
a = 'a';

% [D] Realign to correct motion
r = 'r';

% [E] Coregister and segment

% [G] Debone by removing regions that are not c1, c2 or c3
xDebone_Expr = []; % leave empty for default 'i1+i2+i3>0.5'
b = 'b';

% [U] Reslice tissue probability maps to functional space

% [V] Compute tissue depth map

% [F] Erode c2 and c2 images
xErode_Keep = []; % leave empty for default 0.1, i.e. keep 10% of volume

% [H] [RECOMMENDED] Compute initial DVARS and DVARS spikes
xDvars_do = shiIf(isempty(xDvars_do), true, xDvars_do);
xDvars_SpikeThres = []; % leave empty to use FDR correction, or specify an absolute threshold (e.g. 2)
Dvars = 'Dvars';
SpikeDvars = 'SpikeDvars';

% [I] [RECOMMENDED] Despike voxel by voxel
xDespike_Method = shiIf(isempty(xDespike_Method), '', xDespike_Method); % leave empty for default AFNI '3dDespike', or specify 'ArtRepair'
xDespike_Parameter = {[],[],[]}; % {c1,c2,cOrder} for 3dDespike (default = {2.5, 4.0, nImg/30}), or {WinSize,Thres} for ArtRepair (default = {17, 4}, WinSize must be odd integer)
k = 'k';

% [J] Detrend
xDetrend_Order = []; % leave empty for default round(1+xTr*length(func)/150)
d = 'd';
%%%% prefix Mean will be the same as in [D] realign, HARDWIRED

% [K] Extract global, WM and CSF signals
Nui2 = 'Nui2';
Nui3 = 'Nui3';

% [L] [RECOMMENDED] Interpolate spike timepoints [RECOMMENDED OVER CENSORING. see [P]]
xInterpolate_Method = shiIf(isempty(xInterpolate_Method), '', xInterpolate_Method); % leave empty for default 'Lomb' (Lomb-Scargle periodogram), or enter method for interp1.m ('linear', 'cubic', 'spline', etc)
l = 'l';

% [M] Filter images
% [N] Filter covariates
xFilter_HighCutoff = []; % leave empty for default 0.01
xFilter_LowCutoff = []; % leave empty for default 0.1
xFilter_Nui2_CovFiltInd = []; % leave empty for default filtering all covariates, or 1 * nCov vector of TRUE/FALSE indicating which ones to be filtered
xFilter_Nui3_CovFiltInd = []; % leave empty for default filtering all covariates, or 1 * nCov vector of TRUE/FALSE indicating which ones to be filtered
f = 'f';

% [O] Regress out covariates
v2 = 'v2';
v3 = 'v3';

% [P] [OPTIONAL] Mark spike images for censoring [ONLY IF NOT INTERPOLATED ALREADY. INTERPOLATION RECOMMENDED. see [L]]
Censor = 'Censor';

% [Q] [RECOMMENDED] Smooth
s4 = 's4';
s8 = 's8';

% [R] [RECOMMENDED] Normalize
xNormalize_VoxelSize = shiIf(isempty(xNormalize_VoxelSize), [], xNormalize_VoxelSize); % leave empty for default [3 3 3] (mm)
w = 'w';

% [S] Unwarp atlases
u = 'u';

% [T] Extract time series and compute adjacency matrix
xAdjMat_XtrSummFunc = []; % leave empty for default 'mean', or specify 'med' or 'eig'
xAdjMat_CorrMethod = []; % leave empty for default 'pearson' (same as @(x)atanh(corr(x,'type','pearson')), or specify similarly 'spearman', or any function handle
Adj = 'Adj';

% [W] Compute final DVARS
%%%% all inputs will be the same as in [H] (initial DVARS)

% [X] Preproc summary
PreprocSumm = 'PreprocSumm';


%% ========================================================================

% xxChkFile = fullfile(pth,[anat,nii]);
% assert(exist(xxChkFile,'file'),sprintf('cannot find %s',xxChkFile));
% 
% for xxi = 1:length(func)
%     xxChkFile = fullfile(pth,[func{xxi},nii]);
%     assert(exist(xxChkFile,'file'),sprintf('cannot find %s',xxChkFile));
% end

atlas       = shiStrConcat('atlas_',atlas);
func1       = func{1}     ;
m           = 'm'         ;
c1          = 'c1'        ;
c2          = 'c2'        ;
c3          = 'c3'        ;
rc1         = 'rc1'       ;
rc2         = 'rc2'       ;
rc3         = 'rc3'       ;
ec1         = 'ec1'       ;
ec2         = 'ec2'       ;
ec3         = 'ec3'       ;
iy          = 'iy'        ;
MaskDebone  = 'MaskDebone';
Mean        = 'Mean'      ;
Resliced    = 'Resliced'  ;
Rp          = 'Rp'        ;
Spike       = 'Spike'     ;
Spm         = 'Spm'       ;
TisDep      = 'TisDep'    ;
TisLab      = 'TisLab'    ;
Underline   = '_'         ;
y           = 'y'         ;
BiasField   = 'BiasField' ;
seg8        = 'seg8'      ;
label       = 'label'     ;

nii = '.nii';
txt = '.txt';
mat = '.mat';

ANAT = shiStrConcat( pth, filesep, anat, nii);
FUNC = shiStrConcat( pth, filesep, func, nii);
ATLAS = shiStrConcat( pth, filesep, atlas, nii);
CUSTOMCOV = shiIf( ~isempty(CustomCov), shiStrConcat( pth, filesep, CustomCov, txt ), {} );
CUSTOMSPIKE = shiIf( ~isempty(CustomSpike), shiStrConcat( pth, filesep, CustomSpike, txt ), {} );

manat                      =                             shiStrConcat( pth, filesep, m,anat                                                    , nii );
wanat                      =                             shiStrConcat( pth, filesep, w,anat                                                    , nii );
wmanat                     =                             shiStrConcat( pth, filesep, w,m,anat                                                  , nii );
c1anat                     =                             shiStrConcat( pth, filesep, c1,anat                                                   , nii );
c2anat                     =                             shiStrConcat( pth, filesep, c2,anat                                                   , nii );
c3anat                     =                             shiStrConcat( pth, filesep, c3,anat                                                   , nii );
wc1anat                    =                             shiStrConcat( pth, filesep, w,c1,anat                                                 , nii );
wc2anat                    =                             shiStrConcat( pth, filesep, w,c2,anat                                                 , nii );
wc3anat                    =                             shiStrConcat( pth, filesep, w,c3,anat                                                 , nii );
mwc1anat                   =                             shiStrConcat( pth, filesep, m,w,c1,anat                                               , nii );
mwc2anat                   =                             shiStrConcat( pth, filesep, m,w,c2,anat                                               , nii );
mwc3anat                   =                             shiStrConcat( pth, filesep, m,w,c3,anat                                               , nii );
rc1anat                    =                             shiStrConcat( pth, filesep, rc1,anat                                                  , nii );
rc2anat                    =                             shiStrConcat( pth, filesep, rc2,anat                                                  , nii );
rc3anat                    =                             shiStrConcat( pth, filesep, rc3,anat                                                  , nii );
Resliced_c1anat            =                             shiStrConcat( pth, filesep, Resliced,Underline,c1,anat                                , nii );
Resliced_c2anat            =                             shiStrConcat( pth, filesep, Resliced,Underline,c2,anat                                , nii );
Resliced_c3anat            =                             shiStrConcat( pth, filesep, Resliced,Underline,c3,anat                                , nii );
TisDep_Resliced_c1anat     =                             shiStrConcat( pth, filesep, TisDep,Underline,Resliced,Underline,c1,anat               , nii );
TisLab_Resliced_c1anat     =                             shiStrConcat( pth, filesep, TisLab,Underline,Resliced,Underline,c1,anat               , nii );
ec1_TisDep_Resliced_c1anat =                             shiStrConcat( pth, filesep, ec1,Underline,TisDep,Underline,Resliced,Underline,c1,anat , nii );
ec2_TisDep_Resliced_c1anat =                             shiStrConcat( pth, filesep, ec2,Underline,TisDep,Underline,Resliced,Underline,c1,anat , nii );
ec3_TisDep_Resliced_c1anat =                             shiStrConcat( pth, filesep, ec3,Underline,TisDep,Underline,Resliced,Underline,c1,anat , nii );
y_anat                     =                             shiStrConcat( pth, filesep, y,Underline,anat                                          , nii );
iy_anat                    =                             shiStrConcat( pth, filesep, iy,Underline,anat                                         , nii );
anat_seg8                  =                             shiStrConcat( pth, filesep, anat,Underline,seg8                                       , mat );
BiasField_anat             =                             shiStrConcat( pth, filesep, BiasField,Underline,anat                                  , nii );
uatlas                     =                             shiStrConcat( pth, filesep, u,atlas                                                   , nii );
atlas_label                =                             shiStrConcat( pth, filesep, atlas,Underline,label                                     , txt );
MaskDebone_rafunc1         =                             shiStrConcat( pth, filesep, MaskDebone,Underline,r,a,func1                            , nii );
Resliced_brafunc1          =                             shiStrConcat( pth, filesep, Resliced,Underline,b,r,a,func1                            , nii );
Mean_func1                 =                             shiStrConcat( pth, filesep, Mean,Underline,func1                                      , nii );
Mean_afunc1                =                             shiStrConcat( pth, filesep, Mean,Underline,a,func1                                    , nii );
Mean_brafunc1              =                             shiStrConcat( pth, filesep, Mean,Underline,b,r,a,func1                                , nii );
Mean_kbrafunc1             =                             shiStrConcat( pth, filesep, Mean,Underline,k,b,r,a,func1                              , nii );
Adj_v2dbrafunc1_uatlas     =                             shiStrConcat( pth, filesep, Adj,Underline,v2,d,b,r,a,func1,Underline,u,atlas          , mat );
Adj_v2dkbrafunc1_uatlas    =                             shiStrConcat( pth, filesep, Adj,Underline,v2,d,k,b,r,a,func1,Underline,u,atlas        , mat );
Adj_v2fdbrafunc1_uatlas    =                             shiStrConcat( pth, filesep, Adj,Underline,v2,f,d,b,r,a,func1,Underline,u,atlas        , mat );
Adj_v2fdkbrafunc1_uatlas   =                             shiStrConcat( pth, filesep, Adj,Underline,v2,f,d,k,b,r,a,func1,Underline,u,atlas      , mat );
Adj_v2fldbrafunc1_uatlas   =                             shiStrConcat( pth, filesep, Adj,Underline,v2,f,l,d,b,r,a,func1,Underline,u,atlas      , mat );
Adj_v2fldkbrafunc1_uatlas  =                             shiStrConcat( pth, filesep, Adj,Underline,v2,f,l,d,k,b,r,a,func1,Underline,u,atlas    , mat );
Adj_v2ldbrafunc1_uatlas    =                             shiStrConcat( pth, filesep, Adj,Underline,v2,l,d,b,r,a,func1,Underline,u,atlas        , mat );
Adj_v2ldkbrafunc1_uatlas   =                             shiStrConcat( pth, filesep, Adj,Underline,v2,l,d,k,b,r,a,func1,Underline,u,atlas      , mat );
Adj_v3dbrafunc1_uatlas     =                             shiStrConcat( pth, filesep, Adj,Underline,v3,d,b,r,a,func1,Underline,u,atlas          , mat );
Adj_v3dkbrafunc1_uatlas    =                             shiStrConcat( pth, filesep, Adj,Underline,v3,d,k,b,r,a,func1,Underline,u,atlas        , mat );
Adj_v3fdbrafunc1_uatlas    =                             shiStrConcat( pth, filesep, Adj,Underline,v3,f,d,b,r,a,func1,Underline,u,atlas        , mat );
Adj_v3fdkbrafunc1_uatlas   =                             shiStrConcat( pth, filesep, Adj,Underline,v3,f,d,k,b,r,a,func1,Underline,u,atlas      , mat );
Adj_v3fldbrafunc1_uatlas   =                             shiStrConcat( pth, filesep, Adj,Underline,v3,f,l,d,b,r,a,func1,Underline,u,atlas      , mat );
Adj_v3fldkbrafunc1_uatlas  =                             shiStrConcat( pth, filesep, Adj,Underline,v3,f,l,d,k,b,r,a,func1,Underline,u,atlas    , mat );
Adj_v3ldbrafunc1_uatlas    =                             shiStrConcat( pth, filesep, Adj,Underline,v3,l,d,b,r,a,func1,Underline,u,atlas        , mat );
Adj_v3ldkbrafunc1_uatlas   =                             shiStrConcat( pth, filesep, Adj,Underline,v3,l,d,k,b,r,a,func1,Underline,u,atlas      , mat );
Dvars_brafunc1             = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,b,r,a,func1                               , txt ), {} );
Dvars_v2dbrafunc1          = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,d,b,r,a,func1                          , txt ), {} );
Dvars_v2dkbrafunc1         = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,d,k,b,r,a,func1                        , txt ), {} );
Dvars_v2fdbrafunc1         = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,f,d,b,r,a,func1                        , txt ), {} );
Dvars_v2fdkbrafunc1        = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,f,d,k,b,r,a,func1                      , txt ), {} );
Dvars_v2fldbrafunc1        = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,f,l,d,b,r,a,func1                      , txt ), {} );
Dvars_v2fldkbrafunc1       = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,f,l,d,k,b,r,a,func1                    , txt ), {} );
Dvars_v2ldbrafunc1         = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,l,d,b,r,a,func1                        , txt ), {} );
Dvars_v2ldkbrafunc1        = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v2,l,d,k,b,r,a,func1                      , txt ), {} );
Dvars_v3dbrafunc1          = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,d,b,r,a,func1                          , txt ), {} );
Dvars_v3dkbrafunc1         = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,d,k,b,r,a,func1                        , txt ), {} );
Dvars_v3fdbrafunc1         = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,f,d,b,r,a,func1                        , txt ), {} );
Dvars_v3fdkbrafunc1        = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,f,d,k,b,r,a,func1                      , txt ), {} );
Dvars_v3fldbrafunc1        = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,f,l,d,b,r,a,func1                      , txt ), {} );
Dvars_v3fldkbrafunc1       = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,f,l,d,k,b,r,a,func1                    , txt ), {} );
Dvars_v3ldbrafunc1         = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,l,d,b,r,a,func1                        , txt ), {} );
Dvars_v3ldkbrafunc1        = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Dvars,Underline,v3,l,d,k,b,r,a,func1                      , txt ), {} );
SpikeDvars_brafunc1        = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,b,r,a,func1                         , txt ), {} );
SpikeDvars_v2dbrafunc1     = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,d,b,r,a,func1                    , txt ), {} );
SpikeDvars_v2dkbrafunc1    = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,d,k,b,r,a,func1                  , txt ), {} );
SpikeDvars_v2fdbrafunc1    = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,f,d,b,r,a,func1                  , txt ), {} );
SpikeDvars_v2fdkbrafunc1   = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,f,d,k,b,r,a,func1                , txt ), {} );
SpikeDvars_v2fldbrafunc1   = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,f,l,d,b,r,a,func1                , txt ), {} );
SpikeDvars_v2fldkbrafunc1  = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,f,l,d,k,b,r,a,func1              , txt ), {} );
SpikeDvars_v2ldbrafunc1    = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,l,d,b,r,a,func1                  , txt ), {} );
SpikeDvars_v2ldkbrafunc1   = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v2,l,d,k,b,r,a,func1                , txt ), {} );
SpikeDvars_v3dbrafunc1     = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,d,b,r,a,func1                    , txt ), {} );
SpikeDvars_v3dkbrafunc1    = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,d,k,b,r,a,func1                  , txt ), {} );
SpikeDvars_v3fdbrafunc1    = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,f,d,b,r,a,func1                  , txt ), {} );
SpikeDvars_v3fdkbrafunc1   = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,f,d,k,b,r,a,func1                , txt ), {} );
SpikeDvars_v3fldbrafunc1   = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,f,l,d,b,r,a,func1                , txt ), {} );
SpikeDvars_v3fldkbrafunc1  = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,f,l,d,k,b,r,a,func1              , txt ), {} );
SpikeDvars_v3ldbrafunc1    = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,l,d,b,r,a,func1                  , txt ), {} );
SpikeDvars_v3ldkbrafunc1   = shiIf( xDvars_do,           shiStrConcat( pth, filesep, Spike,Dvars,Underline,v3,l,d,k,b,r,a,func1                , txt ), {} );
Rp_func1                   =                             shiStrConcat( pth, filesep, Rp,Underline,func1                                        , txt );
Rp_afunc1                  =                             shiStrConcat( pth, filesep, Rp,Underline,a,func1                                      , txt );
Fd_Rp_func1                =                             shiStrConcat( pth, filesep, Fd,Underline,Rp,Underline,func1                           , txt );
SpikeFd_Rp_func1           =                             shiStrConcat( pth, filesep, Spike,Fd,Underline,Rp,Underline,func1                     , txt );
Mot24_Rp_func1             =                             shiStrConcat( pth, filesep, Mot24,Underline,Rp,Underline,func1                        , txt );
fMot24_Rp_func1            =                             shiStrConcat( pth, filesep, f,Mot24,Underline,Rp,Underline,func1                      , txt );
Nui2_dbrafunc1             =                             shiStrConcat( pth, filesep, Nui2,Underline,d,b,r,a,func1                              , txt );
Nui2_dkbrafunc1            =                             shiStrConcat( pth, filesep, Nui2,Underline,d,k,b,r,a,func1                            , txt );
Nui3_dbrafunc1             =                             shiStrConcat( pth, filesep, Nui3,Underline,d,b,r,a,func1                              , txt );
Nui3_dkbrafunc1            =                             shiStrConcat( pth, filesep, Nui3,Underline,d,k,b,r,a,func1                            , txt );
fNui2_dbrafunc1            =                             shiStrConcat( pth, filesep, f,Nui2,Underline,d,b,r,a,func1                            , txt );
fNui2_dkbrafunc1           =                             shiStrConcat( pth, filesep, f,Nui2,Underline,d,k,b,r,a,func1                          , txt );
fNui3_dbrafunc1            =                             shiStrConcat( pth, filesep, f,Nui3,Underline,d,b,r,a,func1                            , txt );
fNui3_dkbrafunc1           =                             shiStrConcat( pth, filesep, f,Nui3,Underline,d,k,b,r,a,func1                          , txt );
fCustomCov                 = shiIf( ~isempty(CustomCov), shiStrConcat( pth, filesep, f,CustomCov                                               , txt ), {} );
Censor_brafunc1            =                             shiStrConcat( pth, filesep, Censor,Underline,b,r,a,func1                              , txt );
v2Spm_dbrafunc1            =                             shiStrConcat( pth, filesep, v2,Spm,Underline,d,b,r,a,func1                            , mat );
v2Spm_dkbrafunc1           =                             shiStrConcat( pth, filesep, v2,Spm,Underline,d,k,b,r,a,func1                          , mat );
v2Spm_fdbrafunc1           =                             shiStrConcat( pth, filesep, v2,Spm,Underline,f,d,b,r,a,func1                          , mat );
v2Spm_fdkbrafunc1          =                             shiStrConcat( pth, filesep, v2,Spm,Underline,f,d,k,b,r,a,func1                        , mat );
v2Spm_fldbrafunc1          =                             shiStrConcat( pth, filesep, v2,Spm,Underline,f,l,d,b,r,a,func1                        , mat );
v2Spm_fldkbrafunc1         =                             shiStrConcat( pth, filesep, v2,Spm,Underline,f,l,d,k,b,r,a,func1                      , mat );
v2Spm_ldbrafunc1           =                             shiStrConcat( pth, filesep, v2,Spm,Underline,l,d,b,r,a,func1                          , mat );
v2Spm_ldkbrafunc1          =                             shiStrConcat( pth, filesep, v2,Spm,Underline,l,d,k,b,r,a,func1                        , mat );
v3Spm_dbrafunc1            =                             shiStrConcat( pth, filesep, v3,Spm,Underline,d,b,r,a,func1                            , mat );
v3Spm_dkbrafunc1           =                             shiStrConcat( pth, filesep, v3,Spm,Underline,d,k,b,r,a,func1                          , mat );
v3Spm_fdbrafunc1           =                             shiStrConcat( pth, filesep, v3,Spm,Underline,f,d,b,r,a,func1                          , mat );
v3Spm_fdkbrafunc1          =                             shiStrConcat( pth, filesep, v3,Spm,Underline,f,d,k,b,r,a,func1                        , mat );
v3Spm_fldbrafunc1          =                             shiStrConcat( pth, filesep, v3,Spm,Underline,f,l,d,b,r,a,func1                        , mat );
v3Spm_fldkbrafunc1         =                             shiStrConcat( pth, filesep, v3,Spm,Underline,f,l,d,k,b,r,a,func1                      , mat );
v3Spm_ldbrafunc1           =                             shiStrConcat( pth, filesep, v3,Spm,Underline,l,d,b,r,a,func1                          , mat );
v3Spm_ldkbrafunc1          =                             shiStrConcat( pth, filesep, v3,Spm,Underline,l,d,k,b,r,a,func1                        , mat );
rfunc                      =                             shiStrConcat( pth, filesep, r,func                                                    , nii );
afunc                      =                             shiStrConcat( pth, filesep, a,func                                                    , nii );
rafunc                     =                             shiStrConcat( pth, filesep, r,a,func                                                  , nii );
brafunc                    =                             shiStrConcat( pth, filesep, b,r,a,func                                                , nii );
kbrafunc                   =                             shiStrConcat( pth, filesep, k,b,r,a,func                                              , nii );
dbrafunc                   =                             shiStrConcat( pth, filesep, d,b,r,a,func                                              , nii );
dkbrafunc                  =                             shiStrConcat( pth, filesep, d,k,b,r,a,func                                            , nii );
ldbrafunc                  =                             shiStrConcat( pth, filesep, l,d,b,r,a,func                                            , nii );
ldkbrafunc                 =                             shiStrConcat( pth, filesep, l,d,k,b,r,a,func                                          , nii );
fdbrafunc                  =                             shiStrConcat( pth, filesep, f,d,b,r,a,func                                            , nii );
fdkbrafunc                 =                             shiStrConcat( pth, filesep, f,d,k,b,r,a,func                                          , nii );
fldbrafunc                 =                             shiStrConcat( pth, filesep, f,l,d,b,r,a,func                                          , nii );
fldkbrafunc                =                             shiStrConcat( pth, filesep, f,l,d,k,b,r,a,func                                        , nii );
v2dbrafunc                 =                             shiStrConcat( pth, filesep, v2,d,b,r,a,func                                           , nii );
v2dkbrafunc                =                             shiStrConcat( pth, filesep, v2,d,k,b,r,a,func                                         , nii );
v2fdbrafunc                =                             shiStrConcat( pth, filesep, v2,f,d,b,r,a,func                                         , nii );
v2fdkbrafunc               =                             shiStrConcat( pth, filesep, v2,f,d,k,b,r,a,func                                       , nii );
v2fldbrafunc               =                             shiStrConcat( pth, filesep, v2,f,l,d,b,r,a,func                                       , nii );
v2fldkbrafunc              =                             shiStrConcat( pth, filesep, v2,f,l,d,k,b,r,a,func                                     , nii );
v2ldbrafunc                =                             shiStrConcat( pth, filesep, v2,l,d,b,r,a,func                                         , nii );
v2ldkbrafunc               =                             shiStrConcat( pth, filesep, v2,l,d,k,b,r,a,func                                       , nii );
v3dbrafunc                 =                             shiStrConcat( pth, filesep, v3,d,b,r,a,func                                           , nii );
v3dkbrafunc                =                             shiStrConcat( pth, filesep, v3,d,k,b,r,a,func                                         , nii );
v3fdbrafunc                =                             shiStrConcat( pth, filesep, v3,f,d,b,r,a,func                                         , nii );
v3fdkbrafunc               =                             shiStrConcat( pth, filesep, v3,f,d,k,b,r,a,func                                       , nii );
v3fldbrafunc               =                             shiStrConcat( pth, filesep, v3,f,l,d,b,r,a,func                                       , nii );
v3fldkbrafunc              =                             shiStrConcat( pth, filesep, v3,f,l,d,k,b,r,a,func                                     , nii );
v3ldbrafunc                =                             shiStrConcat( pth, filesep, v3,l,d,b,r,a,func                                         , nii );
v3ldkbrafunc               =                             shiStrConcat( pth, filesep, v3,l,d,k,b,r,a,func                                       , nii );
s4v2dbrafunc               =                             shiStrConcat( pth, filesep, s4,v2,d,b,r,a,func                                        , nii );
s4v2dkbrafunc              =                             shiStrConcat( pth, filesep, s4,v2,d,k,b,r,a,func                                      , nii );
s4v2fdbrafunc              =                             shiStrConcat( pth, filesep, s4,v2,f,d,b,r,a,func                                      , nii );
s4v2fdkbrafunc             =                             shiStrConcat( pth, filesep, s4,v2,f,d,k,b,r,a,func                                    , nii );
s4v2fldbrafunc             =                             shiStrConcat( pth, filesep, s4,v2,f,l,d,b,r,a,func                                    , nii );
s4v2fldkbrafunc            =                             shiStrConcat( pth, filesep, s4,v2,f,l,d,k,b,r,a,func                                  , nii );
s4v2ldbrafunc              =                             shiStrConcat( pth, filesep, s4,v2,l,d,b,r,a,func                                      , nii );
s4v2ldkbrafunc             =                             shiStrConcat( pth, filesep, s4,v2,l,d,k,b,r,a,func                                    , nii );
s4v3dbrafunc               =                             shiStrConcat( pth, filesep, s4,v3,d,b,r,a,func                                        , nii );
s4v3dkbrafunc              =                             shiStrConcat( pth, filesep, s4,v3,d,k,b,r,a,func                                      , nii );
s4v3fdbrafunc              =                             shiStrConcat( pth, filesep, s4,v3,f,d,b,r,a,func                                      , nii );
s4v3fdkbrafunc             =                             shiStrConcat( pth, filesep, s4,v3,f,d,k,b,r,a,func                                    , nii );
s4v3fldbrafunc             =                             shiStrConcat( pth, filesep, s4,v3,f,l,d,b,r,a,func                                    , nii );
s4v3fldkbrafunc            =                             shiStrConcat( pth, filesep, s4,v3,f,l,d,k,b,r,a,func                                  , nii );
s4v3ldbrafunc              =                             shiStrConcat( pth, filesep, s4,v3,l,d,b,r,a,func                                      , nii );
s4v3ldkbrafunc             =                             shiStrConcat( pth, filesep, s4,v3,l,d,k,b,r,a,func                                    , nii );
s8v2dbrafunc               =                             shiStrConcat( pth, filesep, s8,v2,d,b,r,a,func                                        , nii );
s8v2dkbrafunc              =                             shiStrConcat( pth, filesep, s8,v2,d,k,b,r,a,func                                      , nii );
s8v2fdbrafunc              =                             shiStrConcat( pth, filesep, s8,v2,f,d,b,r,a,func                                      , nii );
s8v2fdkbrafunc             =                             shiStrConcat( pth, filesep, s8,v2,f,d,k,b,r,a,func                                    , nii );
s8v2fldbrafunc             =                             shiStrConcat( pth, filesep, s8,v2,f,l,d,b,r,a,func                                    , nii );
s8v2fldkbrafunc            =                             shiStrConcat( pth, filesep, s8,v2,f,l,d,k,b,r,a,func                                  , nii );
s8v2ldbrafunc              =                             shiStrConcat( pth, filesep, s8,v2,l,d,b,r,a,func                                      , nii );
s8v2ldkbrafunc             =                             shiStrConcat( pth, filesep, s8,v2,l,d,k,b,r,a,func                                    , nii );
s8v3dbrafunc               =                             shiStrConcat( pth, filesep, s8,v3,d,b,r,a,func                                        , nii );
s8v3dkbrafunc              =                             shiStrConcat( pth, filesep, s8,v3,d,k,b,r,a,func                                      , nii );
s8v3fdbrafunc              =                             shiStrConcat( pth, filesep, s8,v3,f,d,b,r,a,func                                      , nii );
s8v3fdkbrafunc             =                             shiStrConcat( pth, filesep, s8,v3,f,d,k,b,r,a,func                                    , nii );
s8v3fldbrafunc             =                             shiStrConcat( pth, filesep, s8,v3,f,l,d,b,r,a,func                                    , nii );
s8v3fldkbrafunc            =                             shiStrConcat( pth, filesep, s8,v3,f,l,d,k,b,r,a,func                                  , nii );
s8v3ldbrafunc              =                             shiStrConcat( pth, filesep, s8,v3,l,d,b,r,a,func                                      , nii );
s8v3ldkbrafunc             =                             shiStrConcat( pth, filesep, s8,v3,l,d,k,b,r,a,func                                    , nii );
ws4v2dbrafunc              =                             shiStrConcat( pth, filesep, w,s4,v2,d,b,r,a,func                                      , nii );
ws4v2dkbrafunc             =                             shiStrConcat( pth, filesep, w,s4,v2,d,k,b,r,a,func                                    , nii );
ws4v2fdbrafunc             =                             shiStrConcat( pth, filesep, w,s4,v2,f,d,b,r,a,func                                    , nii );
ws4v2fdkbrafunc            =                             shiStrConcat( pth, filesep, w,s4,v2,f,d,k,b,r,a,func                                  , nii );
ws4v2fldbrafunc            =                             shiStrConcat( pth, filesep, w,s4,v2,f,l,d,b,r,a,func                                  , nii );
ws4v2fldkbrafunc           =                             shiStrConcat( pth, filesep, w,s4,v2,f,l,d,k,b,r,a,func                                , nii );
ws4v2ldbrafunc             =                             shiStrConcat( pth, filesep, w,s4,v2,l,d,b,r,a,func                                    , nii );
ws4v2ldkbrafunc            =                             shiStrConcat( pth, filesep, w,s4,v2,l,d,k,b,r,a,func                                  , nii );
ws4v3dbrafunc              =                             shiStrConcat( pth, filesep, w,s4,v3,d,b,r,a,func                                      , nii );
ws4v3dkbrafunc             =                             shiStrConcat( pth, filesep, w,s4,v3,d,k,b,r,a,func                                    , nii );
ws4v3fdbrafunc             =                             shiStrConcat( pth, filesep, w,s4,v3,f,d,b,r,a,func                                    , nii );
ws4v3fdkbrafunc            =                             shiStrConcat( pth, filesep, w,s4,v3,f,d,k,b,r,a,func                                  , nii );
ws4v3fldbrafunc            =                             shiStrConcat( pth, filesep, w,s4,v3,f,l,d,b,r,a,func                                  , nii );
ws4v3fldkbrafunc           =                             shiStrConcat( pth, filesep, w,s4,v3,f,l,d,k,b,r,a,func                                , nii );
ws4v3ldbrafunc             =                             shiStrConcat( pth, filesep, w,s4,v3,l,d,b,r,a,func                                    , nii );
ws4v3ldkbrafunc            =                             shiStrConcat( pth, filesep, w,s4,v3,l,d,k,b,r,a,func                                  , nii );
ws8v2dbrafunc              =                             shiStrConcat( pth, filesep, w,s8,v2,d,b,r,a,func                                      , nii );
ws8v2dkbrafunc             =                             shiStrConcat( pth, filesep, w,s8,v2,d,k,b,r,a,func                                    , nii );
ws8v2fdbrafunc             =                             shiStrConcat( pth, filesep, w,s8,v2,f,d,b,r,a,func                                    , nii );
ws8v2fdkbrafunc            =                             shiStrConcat( pth, filesep, w,s8,v2,f,d,k,b,r,a,func                                  , nii );
ws8v2fldbrafunc            =                             shiStrConcat( pth, filesep, w,s8,v2,f,l,d,b,r,a,func                                  , nii );
ws8v2fldkbrafunc           =                             shiStrConcat( pth, filesep, w,s8,v2,f,l,d,k,b,r,a,func                                , nii );
ws8v2ldbrafunc             =                             shiStrConcat( pth, filesep, w,s8,v2,l,d,b,r,a,func                                    , nii );
ws8v2ldkbrafunc            =                             shiStrConcat( pth, filesep, w,s8,v2,l,d,k,b,r,a,func                                  , nii );
ws8v3dbrafunc              =                             shiStrConcat( pth, filesep, w,s8,v3,d,b,r,a,func                                      , nii );
ws8v3dkbrafunc             =                             shiStrConcat( pth, filesep, w,s8,v3,d,k,b,r,a,func                                    , nii );
ws8v3fdbrafunc             =                             shiStrConcat( pth, filesep, w,s8,v3,f,d,b,r,a,func                                    , nii );
ws8v3fdkbrafunc            =                             shiStrConcat( pth, filesep, w,s8,v3,f,d,k,b,r,a,func                                  , nii );
ws8v3fldbrafunc            =                             shiStrConcat( pth, filesep, w,s8,v3,f,l,d,b,r,a,func                                  , nii );
ws8v3fldkbrafunc           =                             shiStrConcat( pth, filesep, w,s8,v3,f,l,d,k,b,r,a,func                                , nii );
ws8v3ldbrafunc             =                             shiStrConcat( pth, filesep, w,s8,v3,l,d,b,r,a,func                                    , nii );
ws8v3ldkbrafunc            =                             shiStrConcat( pth, filesep, w,s8,v3,l,d,k,b,r,a,func                                  , nii );
wv2dbrafunc                =                             shiStrConcat( pth, filesep, w,v2,d,b,r,a,func                                         , nii );
wv2dkbrafunc               =                             shiStrConcat( pth, filesep, w,v2,d,k,b,r,a,func                                       , nii );
wv2fdbrafunc               =                             shiStrConcat( pth, filesep, w,v2,f,d,b,r,a,func                                       , nii );
wv2fdkbrafunc              =                             shiStrConcat( pth, filesep, w,v2,f,d,k,b,r,a,func                                     , nii );
wv2fldbrafunc              =                             shiStrConcat( pth, filesep, w,v2,f,l,d,b,r,a,func                                     , nii );
wv2fldkbrafunc             =                             shiStrConcat( pth, filesep, w,v2,f,l,d,k,b,r,a,func                                   , nii );
wv2ldbrafunc               =                             shiStrConcat( pth, filesep, w,v2,l,d,b,r,a,func                                       , nii );
wv2ldkbrafunc              =                             shiStrConcat( pth, filesep, w,v2,l,d,k,b,r,a,func                                     , nii );
wv3dbrafunc                =                             shiStrConcat( pth, filesep, w,v3,d,b,r,a,func                                         , nii );
wv3dkbrafunc               =                             shiStrConcat( pth, filesep, w,v3,d,k,b,r,a,func                                       , nii );
wv3fdbrafunc               =                             shiStrConcat( pth, filesep, w,v3,f,d,b,r,a,func                                       , nii );
wv3fdkbrafunc              =                             shiStrConcat( pth, filesep, w,v3,f,d,k,b,r,a,func                                     , nii );
wv3fldbrafunc              =                             shiStrConcat( pth, filesep, w,v3,f,l,d,b,r,a,func                                     , nii );
wv3fldkbrafunc             =                             shiStrConcat( pth, filesep, w,v3,f,l,d,k,b,r,a,func                                   , nii );
wv3ldbrafunc               =                             shiStrConcat( pth, filesep, w,v3,l,d,b,r,a,func                                       , nii );
wv3ldkbrafunc              =                             shiStrConcat( pth, filesep, w,v3,l,d,k,b,r,a,func                                     , nii );
PreprocSumm_v2fldkbrafunc1 =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,f,l,d,k,b,r,a,func1              , mat );
PreprocSumm_v2fldbrafunc1  =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,f,l,d,b,r,a,func1                , mat );
PreprocSumm_v2fdkbrafunc1  =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,f,d,k,b,r,a,func1                , mat );
PreprocSumm_v2fdbrafunc1   =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,f,d,b,r,a,func1                  , mat );
PreprocSumm_v2ldkbrafunc1  =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,l,d,k,b,r,a,func1                , mat );
PreprocSumm_v2ldbrafunc1   =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,l,d,b,r,a,func1                  , mat );
PreprocSumm_v2dkbrafunc1   =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,d,k,b,r,a,func1                  , mat );
PreprocSumm_v2dbrafunc1    =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v2,d,b,r,a,func1                    , mat );
PreprocSumm_v3fldkbrafunc1 =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,f,l,d,k,b,r,a,func1              , mat );
PreprocSumm_v3fldbrafunc1  =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,f,l,d,b,r,a,func1                , mat );
PreprocSumm_v3fdkbrafunc1  =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,f,d,k,b,r,a,func1                , mat );
PreprocSumm_v3fdbrafunc1   =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,f,d,b,r,a,func1                  , mat );
PreprocSumm_v3ldkbrafunc1  =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,l,d,k,b,r,a,func1                , mat );
PreprocSumm_v3ldbrafunc1   =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,l,d,b,r,a,func1                  , mat );
PreprocSumm_v3dkbrafunc1   =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,d,k,b,r,a,func1                  , mat );
PreprocSumm_v3dbrafunc1    =                             shiStrConcat( pth, filesep, PreprocSumm,Underline,v3,d,b,r,a,func1                    , mat );


%% ========================================================================

xxi = 0;

% [A] * 1           motion estimation               % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           % [A] * 1           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRealign( FUNC, [], 'overwrite' );   OUT{xxi,1} = { Rp_func1; };

% [B] * 2           motion expansion spike          % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           % [B] * 2           
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Motion( Rp_func1, xMotion_FdOption, xMotion_FdSpikeThres, Mot24, Fd, SpikeFd );   OUT{xxi,1} = { Mot24_Rp_func1; Fd_Rp_func1; SpikeFd_Rp_func1; };

% [C] * 3           slice timing                    % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           % [C] * 3           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSliceTime( FUNC, xTr, xTa, xSliceOrder, xRefSlice, a, 'overwrite' );   OUT{xxi,1} = { afunc; };

% [D] * 4           realignment                     % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           % [D] * 4           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRealign( afunc, r, 'overwrite' );   OUT{xxi,1} = { rafunc; Mean_afunc1; };

% [E] * 5           coregistration segmentation     % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           % [E] * 5           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSegment( shiSpmPreprocCoregister( ANAT, Mean_afunc1, [], false, [], 'overwrite' ), true, 'overwrite' );   OUT{xxi,1} = { y_anat; iy_anat; c1anat; c2anat; c3anat; manat; };

% [G] * 6           deboning                        % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           % [G] * 6           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSkullStrip( rafunc, c1anat, c2anat, c3anat, xDebone_Expr, true, b, 'overwrite' );   OUT{xxi,1} = { brafunc; MaskDebone_rafunc1; };

% [U] * 7           reslicing tissue maps           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           % [U] * 7           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmReslice( [brafunc(1), c1anat, c2anat, c3anat], [], [] );   OUT{xxi,1} = { Resliced_c1anat; Resliced_c2anat; Resliced_c3anat; };

% [V] * 8           tissue depth                    % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           % [V] * 8           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocTissueDepth( Resliced_c1anat, Resliced_c2anat, Resliced_c3anat, 'overwrite' );   OUT{xxi,1} = { TisDep_Resliced_c1anat; };

% [F] * 9           tissue eroding                  % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           % [F] * 9           
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocErode( TisDep_Resliced_c1anat, xErode_Keep, 'overwrite' );   OUT{xxi,1} = { ec2_TisDep_Resliced_c1anat; ec3_TisDep_Resliced_c1anat; };

% [H] * 10          dvars1 spike                    % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          % [H] * 10          
xxi = xxi+1;   STEP{xxi,1} = shiIf( ~xDvars_do, @(varargin)disp(''), @() xxStepFunc_Dvars( brafunc, MaskDebone_rafunc1, [], xDvars_SpikeThres, Dvars, SpikeDvars ) );   OUT{xxi,1} = { Dvars_brafunc1; SpikeDvars_brafunc1; };

% [I] * 11          despiking voxelwise             % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          % [I] * 11          
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Despike( xDespike_Method, brafunc, MaskDebone_rafunc1, xDespike_Parameter, k );   OUT{xxi,1} = { kbrafunc; };

% [J] * 12:13       detrending                      % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       % [J] * 12:13       
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocDetrend( kbrafunc, shiIf(~isscalar(xDetrend_Order) || ~isPositiveIntegerValuedNumeric(xDetrend_Order),round(1+xTr*length(func)/150),xDetrend_Order), 0, d, 'overwrite' );   OUT{xxi,1} = { dkbrafunc; Mean_kbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocDetrend(  brafunc, shiIf(~isscalar(xDetrend_Order) || ~isPositiveIntegerValuedNumeric(xDetrend_Order),round(1+xTr*length(func)/150),xDetrend_Order), 0, d, 'overwrite' );   OUT{xxi,1} = {  dbrafunc;  Mean_brafunc1; };

% [K] * 14:17       nuisance extraction             % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       % [K] * 14:17       
xxi = xxi+1;   STEP{xxi,1} = @() writematrix( shiSpmRoiXtr( dkbrafunc, [ ec2_TisDep_Resliced_c1anat; ec3_TisDep_Resliced_c1anat;                     ] ), char(Nui2_dkbrafunc1), 'Delimiter', '\t' );   OUT{xxi,1} = { Nui2_dkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() writematrix( shiSpmRoiXtr(  dbrafunc, [ ec2_TisDep_Resliced_c1anat; ec3_TisDep_Resliced_c1anat;                     ] ), char( Nui2_dbrafunc1), 'Delimiter', '\t' );   OUT{xxi,1} = {  Nui2_dbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() writematrix( shiSpmRoiXtr( dkbrafunc, [ ec2_TisDep_Resliced_c1anat; ec3_TisDep_Resliced_c1anat; MaskDebone_rafunc1; ] ), char(Nui3_dkbrafunc1), 'Delimiter', '\t' );   OUT{xxi,1} = { Nui3_dkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() writematrix( shiSpmRoiXtr(  dbrafunc, [ ec2_TisDep_Resliced_c1anat; ec3_TisDep_Resliced_c1anat; MaskDebone_rafunc1; ] ), char( Nui3_dbrafunc1), 'Delimiter', '\t' );   OUT{xxi,1} = {  Nui3_dbrafunc1; };

% [L] * 18:19       spike interpolation             % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       % [L] * 18:19       
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocInterpolate( dkbrafunc, MaskDebone_rafunc1, [SpikeFd_Rp_func1; SpikeDvars_brafunc1; CUSTOMSPIKE], xTr, xInterpolate_Method, l , 'overwrite' );   OUT{xxi,1} = { ldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocInterpolate(  dbrafunc, MaskDebone_rafunc1, [SpikeFd_Rp_func1; SpikeDvars_brafunc1; CUSTOMSPIKE], xTr, xInterpolate_Method, l , 'overwrite' );   OUT{xxi,1} = {  ldbrafunc; };

% [M] * 20:23       filtering images                % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       % [M] * 20:23       
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocBandPass( ldkbrafunc, xTr, xFilter_HighCutoff, xFilter_LowCutoff, f, 'overwrite' );   OUT{xxi,1} = { fldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocBandPass(  ldbrafunc, xTr, xFilter_HighCutoff, xFilter_LowCutoff, f, 'overwrite' );   OUT{xxi,1} = {  fldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocBandPass(  dkbrafunc, xTr, xFilter_HighCutoff, xFilter_LowCutoff, f, 'overwrite' );   OUT{xxi,1} = {  fdkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocBandPass(   dbrafunc, xTr, xFilter_HighCutoff, xFilter_LowCutoff, f, 'overwrite' );   OUT{xxi,1} = {   fdbrafunc; };

% [N] * 24:27       filtering covariates            % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       % [N] * 24:27       
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_FilterTxt( [Mot24_Rp_func1; Nui2_dkbrafunc1; CustomCov], xTr, xFilter_HighCutoff, xFilter_LowCutoff, xFilter_Nui2_CovFiltInd, f );   OUT{xxi,1} = { fMot24_Rp_func1; fNui2_dkbrafunc1; fCustomCov; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_FilterTxt( [Mot24_Rp_func1;  Nui2_dbrafunc1; CustomCov], xTr, xFilter_HighCutoff, xFilter_LowCutoff, xFilter_Nui2_CovFiltInd, f );   OUT{xxi,1} = { fMot24_Rp_func1;  fNui2_dbrafunc1; fCustomCov; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_FilterTxt( [Mot24_Rp_func1; Nui3_dkbrafunc1; CustomCov], xTr, xFilter_HighCutoff, xFilter_LowCutoff, xFilter_Nui3_CovFiltInd, f );   OUT{xxi,1} = { fMot24_Rp_func1; fNui3_dkbrafunc1; fCustomCov; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_FilterTxt( [Mot24_Rp_func1;  Nui3_dbrafunc1; CustomCov], xTr, xFilter_HighCutoff, xFilter_LowCutoff, xFilter_Nui3_CovFiltInd, f );   OUT{xxi,1} = { fMot24_Rp_func1;  fNui3_dbrafunc1; fCustomCov; };

% [O] * 28:43       regressing out                  % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       % [O] * 28:43       
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut( fldkbrafunc, [], [fMot24_Rp_func1; fNui2_dkbrafunc1; fCustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = { v2fldkbrafunc; v2Spm_fldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(  fldbrafunc, [], [fMot24_Rp_func1;  fNui2_dbrafunc1; fCustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = {  v2fldbrafunc;  v2Spm_fldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(  fdkbrafunc, [], [fMot24_Rp_func1; fNui2_dkbrafunc1; fCustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = {  v2fdkbrafunc;  v2Spm_fdkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(   fdbrafunc, [], [fMot24_Rp_func1;  fNui2_dbrafunc1; fCustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = {   v2fdbrafunc;   v2Spm_fdbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(  ldkbrafunc, [], [ Mot24_Rp_func1;  Nui2_dkbrafunc1;  CustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = {  v2ldkbrafunc;  v2Spm_ldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(   ldbrafunc, [], [ Mot24_Rp_func1;   Nui2_dbrafunc1;  CustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = {   v2ldbrafunc;   v2Spm_ldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(   dkbrafunc, [], [ Mot24_Rp_func1;  Nui2_dkbrafunc1;  CustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = {   v2dkbrafunc;   v2Spm_dkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(    dbrafunc, [], [ Mot24_Rp_func1;   Nui2_dbrafunc1;  CustomCov], [], v2, 'overwrite' );   OUT{xxi,1} = {    v2dbrafunc;    v2Spm_dbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut( fldkbrafunc, [], [fMot24_Rp_func1; fNui3_dkbrafunc1; fCustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = { v3fldkbrafunc; v3Spm_fldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(  fldbrafunc, [], [fMot24_Rp_func1;  fNui3_dbrafunc1; fCustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = {  v3fldbrafunc;  v3Spm_fldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(  fdkbrafunc, [], [fMot24_Rp_func1; fNui3_dkbrafunc1; fCustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = {  v3fdkbrafunc;  v3Spm_fdkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(   fdbrafunc, [], [fMot24_Rp_func1;  fNui3_dbrafunc1; fCustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = {   v3fdbrafunc;   v3Spm_fdbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(  ldkbrafunc, [], [ Mot24_Rp_func1;  Nui3_dkbrafunc1;  CustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = {  v3ldkbrafunc;  v3Spm_ldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(   ldbrafunc, [], [ Mot24_Rp_func1;   Nui3_dbrafunc1;  CustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = {   v3ldbrafunc;   v3Spm_ldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(   dkbrafunc, [], [ Mot24_Rp_func1;  Nui3_dkbrafunc1;  CustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = {   v3dkbrafunc;   v3Spm_dkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocRegressOut(    dbrafunc, [], [ Mot24_Rp_func1;   Nui3_dbrafunc1;  CustomCov], [], v3, 'overwrite' );   OUT{xxi,1} = {    v3dbrafunc;    v3Spm_dbrafunc1; };

% [P] * 44          combining spikes to censors     % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          % [P] * 44          
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Censor( [SpikeFd_Rp_func1; SpikeDvars_brafunc1; CUSTOMSPIKE], Censor_brafunc1 );   OUT{xxi,1} = { Censor_brafunc1 };

% [Q] * 45:76       smoothing                       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       % [Q] * 45:76       
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth( v2fldkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = { s4v2fldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v2fldbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {  s4v2fldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v2fdkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {  s4v2fdkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v2fdbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {   s4v2fdbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v2ldkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {  s4v2ldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v2ldbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {   s4v2ldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v2dkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {   s4v2dkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(    v2dbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {    s4v2dbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth( v3fldkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = { s4v3fldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v3fldbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {  s4v3fldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v3fdkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {  s4v3fdkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v3fdbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {   s4v3fdbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v3ldkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {  s4v3ldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v3ldbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {   s4v3ldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v3dkbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {   s4v3dkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(    v3dbrafunc, [4 4 4], s4, 'overwrite' );   OUT{xxi,1} = {    s4v3dbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth( v2fldkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = { s8v2fldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v2fldbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {  s8v2fldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v2fdkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {  s8v2fdkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v2fdbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {   s8v2fdbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v2ldkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {  s8v2ldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v2ldbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {   s8v2ldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v2dkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {   s8v2dkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(    v2dbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {    s8v2dbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth( v3fldkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = { s8v3fldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v3fldbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {  s8v3fldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v3fdkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {  s8v3fdkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v3fdbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {   s8v3fdbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(  v3ldkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {  s8v3ldkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v3ldbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {   s8v3ldbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(   v3dkbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {   s8v3dkbrafunc; };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocSmooth(    v3dbrafunc, [8 8 8], s8, 'overwrite' );   OUT{xxi,1} = {    s8v3dbrafunc; };

% [R] * 77:125      normalization                   % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      % [R] * 77:125      
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise( s4v2fldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = { ws4v2fldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s4v2fldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws4v2fldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s4v2fdkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws4v2fdkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s4v2fdbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws4v2fdbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s4v2ldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws4v2ldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s4v2ldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws4v2ldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s4v2dkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws4v2dkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    s4v2dbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    ws4v2dbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise( s4v3fldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = { ws4v3fldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s4v3fldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws4v3fldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s4v3fdkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws4v3fdkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s4v3fdbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws4v3fdbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s4v3ldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws4v3ldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s4v3ldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws4v3ldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s4v3dkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws4v3dkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    s4v3dbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    ws4v3dbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise( s8v2fldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = { ws8v2fldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s8v2fldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws8v2fldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s8v2fdkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws8v2fdkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s8v2fdbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws8v2fdbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s8v2ldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws8v2ldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s8v2ldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws8v2ldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s8v2dkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws8v2dkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    s8v2dbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    ws8v2dbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise( s8v3fldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = { ws8v3fldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s8v3fldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws8v3fldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s8v3fdkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws8v3fdkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s8v3fdbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws8v3fdbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(  s8v3ldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {  ws8v3ldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s8v3ldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws8v3ldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   s8v3dkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   ws8v3dkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    s8v3dbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    ws8v3dbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   v2fldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   wv2fldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    v2fldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    wv2fldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    v2fdkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    wv2fdkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(     v2fdbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {     wv2fdbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    v2ldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    wv2ldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(     v2ldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {     wv2ldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(     v2dkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {     wv2dkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(      v2dbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {      wv2dbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(   v3fldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {   wv3fldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    v3fldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    wv3fldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    v3fdkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    wv3fdkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(     v3fdbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {     wv3fdbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    v3ldkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {    wv3ldkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(     v3ldbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {     wv3ldbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(     v3dkbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {     wv3dkbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(      v3dbrafunc, y_anat, xNormalize_VoxelSize, [], [], w, 'overwrite' );   OUT{xxi,1} = {      wv3dbrafunc };
xxi = xxi+1;   STEP{xxi,1} = @() shiSpmPreprocNormalise(    [ANAT;manat], y_anat, [1 1 1], [], [], w, 'overwrite' );   OUT{xxi,1} = { wanat; wmanat; };

% [S] * 126         unwarping atlas                 % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         % [S] * 126         
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AtlasUnwarp( ATLAS, iy_anat, u );   OUT{xxi,1} = { uatlas; };

% [T] * 127:142     adjacency matrix                % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     % [T] * 127:142     
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat( v2fldkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod, Adj_v2fldkbrafunc1_uatlas );   OUT{xxi,1} = { Adj_v2fldkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(  v2fldbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,  Adj_v2fldbrafunc1_uatlas );   OUT{xxi,1} = {  Adj_v2fldbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(  v2fdkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,  Adj_v2fdkbrafunc1_uatlas );   OUT{xxi,1} = {  Adj_v2fdkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(   v2fdbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,   Adj_v2fdbrafunc1_uatlas );   OUT{xxi,1} = {   Adj_v2fdbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(  v2ldkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,  Adj_v2ldkbrafunc1_uatlas );   OUT{xxi,1} = {  Adj_v2ldkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(   v2ldbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,   Adj_v2ldbrafunc1_uatlas );   OUT{xxi,1} = {   Adj_v2ldbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(   v2dkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,   Adj_v2dkbrafunc1_uatlas );   OUT{xxi,1} = {   Adj_v2dkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(    v2dbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,    Adj_v2dbrafunc1_uatlas );   OUT{xxi,1} = {    Adj_v2dbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat( v3fldkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod, Adj_v3fldkbrafunc1_uatlas );   OUT{xxi,1} = { Adj_v3fldkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(  v3fldbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,  Adj_v3fldbrafunc1_uatlas );   OUT{xxi,1} = {  Adj_v3fldbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(  v3fdkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,  Adj_v3fdkbrafunc1_uatlas );   OUT{xxi,1} = {  Adj_v3fdkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(   v3fdbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,   Adj_v3fdbrafunc1_uatlas );   OUT{xxi,1} = {   Adj_v3fdbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(  v3ldkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,  Adj_v3ldkbrafunc1_uatlas );   OUT{xxi,1} = {  Adj_v3ldkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(   v3ldbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,   Adj_v3ldbrafunc1_uatlas );   OUT{xxi,1} = {   Adj_v3ldbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(   v3dkbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,   Adj_v3dkbrafunc1_uatlas );   OUT{xxi,1} = {   Adj_v3dkbrafunc1_uatlas };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_AdjMat(    v3dbrafunc, uatlas, atlas_label, Censor_brafunc1, xAdjMat_XtrSummFunc, xAdjMat_CorrMethod,    Adj_v3dbrafunc1_uatlas );   OUT{xxi,1} = {    Adj_v3dbrafunc1_uatlas };

% [W] * 143:158     dvars2                          % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     % [W] * 143:158     
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars( v2fldkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = { Dvars_v2fldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(  v2fldbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {  Dvars_v2fldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(  v2fdkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {  Dvars_v2fdkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(   v2fdbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {   Dvars_v2fdbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(  v2ldkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {  Dvars_v2ldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(   v2ldbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {   Dvars_v2ldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(   v2dkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {   Dvars_v2dkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(    v2dbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {    Dvars_v2dbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars( v3fldkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = { Dvars_v3fldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(  v3fldbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {  Dvars_v3fldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(  v3fdkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {  Dvars_v3fdkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(   v3fdbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {   Dvars_v3fdbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(  v3ldkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {  Dvars_v3ldkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(   v3ldbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {   Dvars_v3ldbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(   v3dkbrafunc, MaskDebone_rafunc1, Mean_kbrafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {   Dvars_v3dkbrafunc1; };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_Dvars(    v3dbrafunc, MaskDebone_rafunc1,  Mean_brafunc1, xDvars_SpikeThres, Dvars, SpikeDvars );   OUT{xxi,1} = {    Dvars_v3dbrafunc1; };

% [X] * 159:174     summary                         % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     % [X] * 159:174     
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc, v2fldkbrafunc, Fd_Rp_func1, Dvars_brafunc1, Dvars_v2fldkbrafunc1, CustomCov, v2Spm_fldkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = { PreprocSumm_v2fldkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,  v2fldbrafunc, Fd_Rp_func1, Dvars_brafunc1,  Dvars_v2fldbrafunc1, CustomCov,  v2Spm_fldbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {  PreprocSumm_v2fldbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc,  v2fdkbrafunc, Fd_Rp_func1, Dvars_brafunc1,  Dvars_v2fdkbrafunc1, CustomCov,  v2Spm_fdkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {  PreprocSumm_v2fdkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,   v2fdbrafunc, Fd_Rp_func1, Dvars_brafunc1,   Dvars_v2fdbrafunc1, CustomCov,   v2Spm_fdbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {   PreprocSumm_v2fdbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc,  v2ldkbrafunc, Fd_Rp_func1, Dvars_brafunc1,  Dvars_v2ldkbrafunc1, CustomCov,  v2Spm_ldkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {  PreprocSumm_v2ldkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,   v2ldbrafunc, Fd_Rp_func1, Dvars_brafunc1,   Dvars_v2ldbrafunc1, CustomCov,   v2Spm_ldbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {   PreprocSumm_v2ldbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc,   v2dkbrafunc, Fd_Rp_func1, Dvars_brafunc1,   Dvars_v2dkbrafunc1, CustomCov,   v2Spm_dkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {   PreprocSumm_v2dkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,    v2dbrafunc, Fd_Rp_func1, Dvars_brafunc1,    Dvars_v2dbrafunc1, CustomCov,    v2Spm_dbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {    PreprocSumm_v2dbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc, v3fldkbrafunc, Fd_Rp_func1, Dvars_brafunc1, Dvars_v3fldkbrafunc1, CustomCov, v3Spm_fldkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = { PreprocSumm_v3fldkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,  v3fldbrafunc, Fd_Rp_func1, Dvars_brafunc1,  Dvars_v3fldbrafunc1, CustomCov,  v3Spm_fldbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {  PreprocSumm_v3fldbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc,  v3fdkbrafunc, Fd_Rp_func1, Dvars_brafunc1,  Dvars_v3fdkbrafunc1, CustomCov,  v3Spm_fdkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {  PreprocSumm_v3fdkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,   v3fdbrafunc, Fd_Rp_func1, Dvars_brafunc1,   Dvars_v3fdbrafunc1, CustomCov,   v3Spm_fdbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {   PreprocSumm_v3fdbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc,  v3ldkbrafunc, Fd_Rp_func1, Dvars_brafunc1,  Dvars_v3ldkbrafunc1, CustomCov,  v3Spm_ldkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {  PreprocSumm_v3ldkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,   v3ldbrafunc, Fd_Rp_func1, Dvars_brafunc1,   Dvars_v3ldbrafunc1, CustomCov,   v3Spm_ldbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {   PreprocSumm_v3ldbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat, dkbrafunc,   v3dkbrafunc, Fd_Rp_func1, Dvars_brafunc1,   Dvars_v3dkbrafunc1, CustomCov,   v3Spm_dkbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {   PreprocSumm_v3dkbrafunc1 };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFunc_PreprocSumm( TisDep_Resliced_c1anat,  dbrafunc,    v3dbrafunc, Fd_Rp_func1, Dvars_brafunc1,    Dvars_v3dbrafunc1, CustomCov,    v3Spm_dbrafunc1, SpikeFd_Rp_func1, SpikeDvars_brafunc1, CUSTOMSPIKE, Censor_brafunc1, PreprocSumm );   OUT{xxi,1} = {    PreprocSumm_v3dbrafunc1 };

% [Y] * 175:177     tar                             % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     % [Y] * 175:177     
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFun_Compress({
    [
    ATLAS; uatlas; atlas_label;
    ];
    [
    ANAT; manat; wanat; wmanat;
    c1anat; c2anat; c3anat; rc1anat; rc2anat; rc3anat; Resliced_c1anat; Resliced_c2anat; Resliced_c3anat;
    wc1anat; wc2anat; wc3anat; mwc1anat; mwc2anat; mwc3anat;
    TisLab_Resliced_c1anat; TisDep_Resliced_c1anat;
    ec1_TisDep_Resliced_c1anat; ec2_TisDep_Resliced_c1anat; ec3_TisDep_Resliced_c1anat;
    y_anat; iy_anat; BiasField_anat; anat_seg8;
    ];
    [
    Mean_func1; Mean_afunc1; Mean_brafunc1; Mean_kbrafunc1; Resliced_brafunc1;
    ];
    },{
    fullfile(pth,'atlas_etc.gz');
    fullfile(pth,'anat_etc.gz');
    fullfile(pth,'other_images.gz');
    });   OUT{xxi,1} = {  };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFun_Compress({
    [
    Dvars_brafunc1      ;
    Dvars_v2dbrafunc1   ;
    Dvars_v2dkbrafunc1  ;
    Dvars_v2fdbrafunc1  ;
    Dvars_v2fdkbrafunc1 ;
    Dvars_v2fldbrafunc1 ;
    Dvars_v2fldkbrafunc1;
    Dvars_v2ldbrafunc1  ;
    Dvars_v2ldkbrafunc1 ;
    Dvars_v3dbrafunc1   ;
    Dvars_v3dkbrafunc1  ;
    Dvars_v3fdbrafunc1  ;
    Dvars_v3fdkbrafunc1 ;
    Dvars_v3fldbrafunc1 ;
    Dvars_v3fldkbrafunc1;
    Dvars_v3ldbrafunc1  ;
    Dvars_v3ldkbrafunc1 ;
    Rp_func1            ;
    Rp_afunc1           ;
    Fd_Rp_func1         ;
    Mot24_Rp_func1      ;
    fMot24_Rp_func1     ;
    Nui2_dbrafunc1      ;
    Nui2_dkbrafunc1     ;
    Nui3_dbrafunc1      ;
    Nui3_dkbrafunc1     ;
    fNui2_dbrafunc1     ;
    fNui2_dkbrafunc1    ;
    fNui3_dbrafunc1     ;
    fNui3_dkbrafunc1    ;
    CUSTOMCOV           ;
    fCustomCov          ;
    v2Spm_dbrafunc1     ;
    v2Spm_dkbrafunc1    ;
    v2Spm_fdbrafunc1    ;
    v2Spm_fdkbrafunc1   ;
    v2Spm_fldbrafunc1   ;
    v2Spm_fldkbrafunc1  ;
    v2Spm_ldbrafunc1    ;
    v2Spm_ldkbrafunc1   ;
    v3Spm_dbrafunc1     ;
    v3Spm_dkbrafunc1    ;
    v3Spm_fdbrafunc1    ;
    v3Spm_fdkbrafunc1   ;
    v3Spm_fldbrafunc1   ;
    v3Spm_fldkbrafunc1  ;
    v3Spm_ldbrafunc1    ;
    v3Spm_ldkbrafunc1   ;
    ];
    [
    SpikeDvars_brafunc1      ;
    SpikeDvars_v2dbrafunc1   ;
    SpikeDvars_v2dkbrafunc1  ;
    SpikeDvars_v2fdbrafunc1  ;
    SpikeDvars_v2fdkbrafunc1 ;
    SpikeDvars_v2fldbrafunc1 ;
    SpikeDvars_v2fldkbrafunc1;
    SpikeDvars_v2ldbrafunc1  ;
    SpikeDvars_v2ldkbrafunc1 ;
    SpikeDvars_v3dbrafunc1   ;
    SpikeDvars_v3dkbrafunc1  ;
    SpikeDvars_v3fdbrafunc1  ;
    SpikeDvars_v3fdkbrafunc1 ;
    SpikeDvars_v3fldbrafunc1 ;
    SpikeDvars_v3fldkbrafunc1;
    SpikeDvars_v3ldbrafunc1  ;
    SpikeDvars_v3ldkbrafunc1 ;
    SpikeFd_Rp_func1         ;
    Censor_brafunc1          ;
    ];
    },{
    fullfile(pth,'Covariates.gz');
    fullfile(pth,'Spikes.gz');
    });   OUT{xxi,1} = {  };
xxi = xxi+1;   STEP{xxi,1} = @() xxStepFun_Compress({
    FUNC               ;
    afunc              ;
    brafunc            ;
    dbrafunc           ;
    dkbrafunc          ;
    fdbrafunc          ;
    fdkbrafunc         ;
    fldbrafunc         ;
    fldkbrafunc        ;
    kbrafunc           ;
    ldbrafunc          ;
    ldkbrafunc         ;
    rafunc             ;
    rfunc              ;
    s4v2dbrafunc       ;
    s4v2dkbrafunc      ;
    s4v2fdbrafunc      ;
    s4v2fdkbrafunc     ;
    s4v2fldbrafunc     ;
    s4v2fldkbrafunc    ;
    s4v2ldbrafunc      ;
    s4v2ldkbrafunc     ;
    s4v3dbrafunc       ;
    s4v3dkbrafunc      ;
    s4v3fdbrafunc      ;
    s4v3fdkbrafunc     ;
    s4v3fldbrafunc     ;
    s4v3fldkbrafunc    ;
    s4v3ldbrafunc      ;
    s4v3ldkbrafunc     ;
    s8v2dbrafunc       ;
    s8v2dkbrafunc      ;
    s8v2fdbrafunc      ;
    s8v2fdkbrafunc     ;
    s8v2fldbrafunc     ;
    s8v2fldkbrafunc    ;
    s8v2ldbrafunc      ;
    s8v2ldkbrafunc     ;
    s8v3dbrafunc       ;
    s8v3dkbrafunc      ;
    s8v3fdbrafunc      ;
    s8v3fdkbrafunc     ;
    s8v3fldbrafunc     ;
    s8v3fldkbrafunc    ;
    s8v3ldbrafunc      ;
    s8v3ldkbrafunc     ;
    v2dbrafunc         ;
    v2dkbrafunc        ;
    v2fdbrafunc        ;
    v2fdkbrafunc       ;
    v2fldbrafunc       ;
    v2fldkbrafunc      ;
    v2ldbrafunc        ;
    v2ldkbrafunc       ;
    v3dbrafunc         ;
    v3dkbrafunc        ;
    v3fdbrafunc        ;
    v3fdkbrafunc       ;
    v3fldbrafunc       ;
    v3fldkbrafunc      ;
    v3ldbrafunc        ;
    v3ldkbrafunc       ;
    ws4v2dbrafunc      ;
    ws4v2dkbrafunc     ;
    ws4v2fdbrafunc     ;
    ws4v2fdkbrafunc    ;
    ws4v2fldbrafunc    ;
    ws4v2fldkbrafunc   ;
    ws4v2ldbrafunc     ;
    ws4v2ldkbrafunc    ;
    ws4v3dbrafunc      ;
    ws4v3dkbrafunc     ;
    ws4v3fdbrafunc     ;
    ws4v3fdkbrafunc    ;
    ws4v3fldbrafunc    ;
    ws4v3fldkbrafunc   ;
    ws4v3ldbrafunc     ;
    ws4v3ldkbrafunc    ;
    ws8v2dbrafunc      ;
    ws8v2dkbrafunc     ;
    ws8v2fdbrafunc     ;
    ws8v2fdkbrafunc    ;
    ws8v2fldbrafunc    ;
    ws8v2fldkbrafunc   ;
    ws8v2ldbrafunc     ;
    ws8v2ldkbrafunc    ;
    ws8v3dbrafunc      ;
    ws8v3dkbrafunc     ;
    ws8v3fdbrafunc     ;
    ws8v3fdkbrafunc    ;
    ws8v3fldbrafunc    ;
    ws8v3fldkbrafunc   ;
    ws8v3ldbrafunc     ;
    ws8v3ldkbrafunc    ;
    wv2dbrafunc        ;
    wv2dkbrafunc       ;
    wv2fdbrafunc       ;
    wv2fdkbrafunc      ;
    wv2fldbrafunc      ;
    wv2fldkbrafunc     ;
    wv2ldbrafunc       ;
    wv2ldkbrafunc      ;
    wv3dbrafunc        ;
    wv3dkbrafunc       ;
    wv3fdbrafunc       ;
    wv3fdkbrafunc      ;
    wv3fldbrafunc      ;
    wv3fldkbrafunc     ;
    wv3ldbrafunc       ;
    wv3ldkbrafunc      ;
    },[]); OUT{xxi,1} = {  }; %#ok<*NASGU> 

%% ========================================================================
%% ========================================================================
%% ========================================================================
%%
for i = xStepInd
    shiDisp({ ...
        shiTime(0), ...
        '', ...
        'Path:', ...
        ['    ',pth], ...
        '', ...
        sprintf('Step %d (of %d):',i,xxi), ...
        ['    ',char(STEP{i})] ...
        });
    STEP{i}();
end

%% ========================================================================
%% ========================================================================
%% ========================================================================

%%
function xxStepFunc_Motion( Rp, FdOption, FdSpikeThres, pfxMot24, pfxFd, pfxSpikeFd )

Rp = char(Rp);
if ~exist('FdOption','var') || isempty(FdOption)
    FdOption = 'Fd_Jenkinson';
end
if ~exist('FdSpikeThres','var') || isempty(FdSpikeThres)
    FdSpikeThres = 0.5;
end

Mot_dat = shiSpmPreprocMotCalc(Rp);
Mot24 = [Mot_dat.Abs, Mot_dat.Fd, Mot_dat.AbsSq, Mot_dat.FdSq];
Fd = Mot_dat.(FdOption);
Spk = Fd > FdSpikeThres;

[pth,nme] = fileparts(Rp);
writematrix( Mot24, fullfile(pth, [ pfxMot24   ,'_',nme, '.txt' ]), 'Delimiter', '\t' );
writematrix( Fd   , fullfile(pth, [ pfxFd      ,'_',nme, '.txt' ]), 'Delimiter', '\t' );
writematrix( Spk  , fullfile(pth, [ pfxSpikeFd ,'_',nme, '.txt' ]), 'Delimiter', '\t' );


%%
function xxStepFunc_Dvars( Img, MaskImg, MeanImg, SpikeThres, pfxDvars, pfxSpikeDvars )

[Dvars,~,SpkDef] = shiSpmPreprocDvarsCalc(Img,MaskImg,MeanImg);

if ~exist('SpikeThres','var') || isempty(SpikeThres) || ( ischar(SpikeThres) && shiStrIncl(lower(SpikeThres),{'def'}) )
    Spk = SpkDef; % p_fwe<0.05, delta%D-vars>5
elseif isnumeric(SpikeThres) && SpikeThres > 0
    Spk = Dvars > SpikeThres;
else
    error('spike threshold must be ''default'', a positive z-threshod, or left empty for ''default''');
end

[pth,nme] = fileparts(Img{1});
writematrix( Dvars, fullfile(pth, [ pfxDvars      ,'_',nme, '.txt' ]), 'delimiter', '\t' );
writematrix( Spk  , fullfile(pth, [ pfxSpikeDvars ,'_',nme, '.txt' ]), 'delimiter', '\t' );


%%
function xxStepFunc_Despike( Method, Img, MaskImg, Param, k )

if matches(lower(Method),'artrepair')
    shiSpmPreprocArtDespike( Img, Param{1}, Param{2}, k, 'overwrite' );
else
    shiSpmPreproc3dDespike( Img, MaskImg, Param{1}, Param{2}, Param{3}, k, 'overwrite' ) ;
end


%%
function xxStepFunc_FilterTxt( CovFile, Tr, HighCutoff, LowCutoff, FiltInd, pfx )

if isempty(CovFile)
    return;
end

[X,indTxt] = deal([]);
for i = 1:numel(CovFile)
    if isempty(CovFile{i})
        continue;
    end
    xx = readmatrix(CovFile{i});
    if isempty(xx)
        continue;
    end
    X = cat(2,X,xx);
    indTxt = cat(2,indTxt,i.*ones(1,size(xx,2)));
end

X_filtered = X;

if ~exist('FiltInd','var') || isempty(FiltInd) || isequal(FiltInd,1) || isequal(FiltInd,true)
    FiltInd = true(size(indTxt));
end
assert(isequal(numel(indTxt),numel(FiltInd)));

for i = 1:size(X,2)
    if FiltInd(i)
        X_filtered(:,i) = shiStatBandPass(X(:,i),1/Tr,HighCutoff,LowCutoff);
    end
end

for i = 1:numel(CovFile)
    if isempty(CovFile{i}) || ~any(indTxt==i)
        continue;
    end
    [pth,nme,ext] = fileparts(CovFile{i});
    writematrix(X_filtered(:,indTxt==i),fullfile(pth,[pfx,nme,ext]),'Delimiter','\t');
end


%%
function xxStepFunc_Censor( SpikeFile, OutputFile )

X = [];

for i = 1:numel(SpikeFile)
    if isempty(SpikeFile{i})
        continue;
    end
    X = cat(2,X,logical(readmatrix(SpikeFile{i})));
end

writematrix(X,char(OutputFile),'Delimiter','\t');


%%
function xxStepFunc_AtlasUnwarp( ATLAS, iy, pfx )

for i = 1:length(ATLAS)

    [pth,xatlas] = fileparts(ATLAS{i});
    xatlas = xatlas(7:end);

    [xxnii,xxtxt] = shiSpmAnatLabel(xatlas);
    copyfile(xxnii,fullfile(pth,['atlas_',xatlas,'.nii']));
    copyfile(xxtxt,fullfile(pth,['atlas_',xatlas,'_label.txt']));
end

shiSpmPreprocNormalise( ATLAS, iy, [1 1 1], 'Deform', 0, pfx, 'overwrite' );


%%
function xxStepFunc_AdjMat( Img, AtlasAll, AtlasLabAll, CensorFile, XtrSummFunc, CorrMethod, MatName )

if isempty(CorrMethod) || (ischar(CorrMethod) && isequal(lower(CorrMethod),'pearson'))
    CorrMethod = @(x)atanh(corr(x,'type','pearson'));
elseif ischar(CorrMethod) && isequal(lower(CorrMethod),'spearman')
    CorrMethod = @(x)atanh(corr(x,'type','spearman'));
end

assert(isa(CorrMethod,'function_handle'),'unknown correlation function');
assert(numel(AtlasAll) == numel(AtlasLabAll),'unequal length of atlas files and atlas text files');
assert(numel(AtlasAll) == numel(MatName),'unequal length of atlas files and output mat files');
CensorFile = char(CensorFile);

if ~isempty(CensorFile)
    Censor = readmatrix(CensorFile);
    if isempty(Censor)
        Censor = false(size(Img));
    else
        Censor = any(Censor,2);
    end
else
    Censor = false(size(Img));
end

allCensor = all(Censor);
noCensor = all(~Censor);
if allCensor
    warning('cannot censor all timepints.');
end

for i = 1:length(AtlasAll)
    X = shiSpmRoiXtrMultilabel(Img,AtlasAll{i},[],AtlasLabAll{i},XtrSummFunc);
    Adj = CorrMethod(X);
    if noCensor
        Adj_censored = Adj;
    elseif allCensor
        Adj_censored = nan(size(Adj));
    else
        Adj_censored = CorrMethod(X(~Censor,:));
    end
    Atlas = AtlasAll{i};
    AtlasLabFile = AtlasLabAll{i};
    AtlasLab = readtable(AtlasLabFile,'NumHeaderLines',0,'Delimiter','\t');
    AtlasLab.Properties.VariableNames = {'Value','Name'};
    TIME = shiTime(0);
    cwd = pwd;
    save(MatName{i},'Adj','Adj_censored','X','Censor','Img','Atlas','AtlasLab','AtlasLabFile','CensorFile','XtrSummFunc','CorrMethod','TIME','cwd');
end


%%
function xxStepFunc_PreprocSumm( ImgTisFile, ImgMinimalFile, ImgFinalFile, FdFile, DvarsFile, DvarsFinalFile, CovCustFile, SpmFile, SpkFdFile, SpkDvarsFile, SpkCustFile, CensFile, pfx )

PercentVoxelSampled = 10; % e.g. 10, for 10 percent of voxels sampled

nVol = length(ImgMinimalFile);

[pth,nme] = fileparts(ImgFinalFile{1});
OUT_NAME = [pfx,'_',nme];
OUT_NAME_MAT = fullfile(pth,[OUT_NAME,'.mat']);
OUT_NAME_FIG = fullfile(pth,[OUT_NAME,'.jpg']);

Fd = abs(xxreadtxt(FdFile,nVol,NaN));
SpkFd = xxreadtxt(SpkFdFile,nVol,false);

Dvars = xxreadtxt(DvarsFile,nVol,NaN);
DvarsFinal = xxreadtxt(DvarsFinalFile,nVol,NaN);
SpkDvars = xxreadtxt(SpkDvarsFile,nVol,false);

CovCust = xxreadtxt(CovCustFile,nVol,NaN);
SpkCust = xxreadtxt(SpkCustFile,nVol,false);

Spm = load(char(SpmFile));

Cens = xxreadtxt(CensFile,nVol,NaN);

nSpikeFd = sum(SpkFd);
nSpikeDvars = sum(SpkDvars);
nSpikeCustom = sum(SpkCust);
nSpikeALL = sum(any([SpkFd,SpkDvars,SpkCust],2));

FdMean = mean(Fd);
FdMax = max(Fd);

AllCov = Spm.SPM.xX.X;
LostDF_Cov = rank(AllCov) - 1;
LostDF_Cens = sum(any(Cens,2));

corrFdDvars1 = corr(Dvars,Fd);
corrFdDvars2 = corr(DvarsFinal,Fd);

ImgMinimal = xxreadnii(ImgMinimalFile,PercentVoxelSampled);
ImgFinal = xxreadnii(ImgFinalFile,PercentVoxelSampled);
ImgTis = xxreadnii(ImgTisFile,PercentVoxelSampled);

idx = ImgTis > 1;
ImgMinimal = ImgMinimal(:,idx);
ImgFinal = ImgFinal(:,idx);
ImgTis = ImgTis(:,idx);

[~,idx] = sort(ImgTis);
ImgMinimal = ImgMinimal(:,idx);
ImgFinal = ImgFinal(:,idx);
ImgTis = ImgTis(:,idx);

ImgMinimal = zscore(ImgMinimal);
ImgFinal = zscore(ImgFinal);

f = figure('visible','off');
colormap spring;

subplot(10,10,1:10);
imagesc([1:nVol,1:nVol]',[ones(nVol,1)*min(Dvars);ones(nVol,1)*max(Dvars)],[Dvars,Dvars]'); set(gca,'YDir','normal'); hold on;
plot((1:nVol)',Dvars,'k-',shiIf(~any(SpkDvars),NaN,find(SpkDvars)),max(Dvars),'bd'); ylim([min(Dvars),max(Dvars)]); xticklabels([]);
ylabel('DVARS');
title(OUT_NAME,' ','interpreter','none'); hold off;

subplot(10,10,11:20);
imagesc([1:nVol,1:nVol]',[ones(nVol,1)*min(Fd);ones(nVol,1)*max(Fd)],[Fd,Fd]'); set(gca,'YDir','normal'); hold on;
plot((1:nVol)',Fd,'k-',shiIf(~any(SpkFd),NaN,find(SpkFd)),max(Fd),'bd'); ylim([min(Fd),max(Fd)]); xticklabels([]);
ylabel('FD');
hold off;

subplot(10,10,21:60);
imagesc((1:nVol)',ImgTis,ImgMinimal'); hold on;
plot([1;nVol],[2,2],'k-',[1;nVol],[3,3],'k-');
yticks(1.5:1:3.5); yticklabels({'GM','WM','CSF'}); ytickangle(90);
xticklabels([]); ylabel('MINIMAL');
hold off;

subplot(10,10,61:100);
imagesc((1:nVol)',ImgTis,ImgFinal'); hold on;
plot([1;nVol],[2,2],'k-',[1;nVol],[3,3],'k-');
yticks(1.5:1:3.5); yticklabels({'GM','WM','CSF'}); ytickangle(90);
ylabel('FULL');
hold off;

set(f, 'PaperUnits', 'inches');
set(f, 'PaperPosition', [0 0 10 7]);
saveas(f, OUT_NAME_FIG);

Input = struct( ...
    'ImgTis'    , {ImgTisFile    }, ...
    'ImgMinimal', {ImgMinimalFile}, ...
    'ImgFinal'  , {ImgFinalFile  }, ...
    'Fd'        , {FdFile        }, ...
    'Dvars'     , {DvarsFile     }, ...
    'DvarsFinal', {DvarsFinalFile}, ...
    'CovCust'   , {CovCustFile   }, ...
    'Spm'       , {SpmFile       }, ...
    'SpkFd'     , {SpkFdFile     }, ...
    'SpkDvars'  , {SpkDvarsFile  }, ...
    'SpkCust'   , {SpkCustFile   }, ...
    'Cens'      , {CensFile      }, ...
    'pfx'       , {pfx           }, ...
    'PercentVoxelSampled', {PercentVoxelSampled} ...
    );
Data = struct( ...
    'Fd'        , {Fd        }, ...
    'SpkFd'     , {SpkFd     }, ...
    'Dvars'     , {Dvars     }, ...
    'SpkDvars'  , {SpkDvars  }, ...
    'DvarsFinal', {DvarsFinal}, ...
    'CovCust'   , {CovCust   }, ...
    'SpkCust'   , {SpkCust   }, ...
    'Cens'      , {Cens      }, ...
    'AllCov'    , {AllCov    } ...
    );
Result = struct( ...
    'nVol'        , {nVol        }, ...
    'nSpikeFd'    , {nSpikeFd    }, ...
    'nSpikeDvars' , {nSpikeDvars }, ...
    'nSpikeCustom', {nSpikeCustom}, ...
    'nSpikeALL'   , {nSpikeALL   }, ...
    'FdMean'      , {FdMean      }, ...
    'FdMax'       , {FdMax       }, ...
    'LostDF_Cov'  , {LostDF_Cov  }, ...
    'LostDF_Cens' , {LostDF_Cens }, ...
    'corrFdDvars1', {corrFdDvars1}, ...
    'corrFdDvars2', {corrFdDvars2}, ...
    'Figure_filename', {OUT_NAME_FIG} ...
    );

save(OUT_NAME_MAT,'Input','Data','Result');


%%
function xxStepFun_Compress( F, gzName )

if isempty(F)
    return;
end

if ~exist('gzName','var') || isempty(gzName)
    for i = 1:length(F)
        if isempty(F{i})
            gzName{i} = '';
        else
            gzName{i} = [F{i}{1},'.gz'];
        end
    end
end

assert(numel(gzName)==numel(F))

for i = 1:length(F)

    fprintf('(%3d/%d) - gzipping %s ...\n',i,length(F),F{i}{1});

    if isempty(F{i})
        continue;
    elseif ~all(cellfun(@(x)exist(x,'file'),F{i}))
        warning('incomplete file list %s', F{i}{1});
        continue;
    end
    if exist(gzName{i},'file')
        warning('%s already exists',gzName{i});
        continue;
    end

    xximg16(F{i});
    tar(gzName{i},F{i});
    if exist(gzName{i},'file')
        cellfun(@delete,F{i})
    end

end


%%
function x = xxreadtxt(x,n,def)
try x = readmatrix(char(x));
catch, x = repmat(def,n,1);
end

function x = xxreadnii(x,PercentVoxelSampled)
x = spm_read_vols(spm_vol(char(x)));
x = reshape(x,prod(size(x,1:3)),size(x,4))';
k = round(100/PercentVoxelSampled);
x = x(:,1:k:end);

function xximg16(f)
for i = 1:length(f)
    [pt,~,xt] = fileparts(f{i});
    if ~matches(lower(xt),{'.nii','.img'}), continue; end
    V = spm_vol(f{i});
    if V.dt(1)~=64, continue; end
    Y = single(spm_read_vols(V));
    V.dt(1) = 16;
    tmp = fullfile(pt,['_tmp',xt]);
    movefile(f{i},tmp);
    spm_write_vol(V,Y);
    assert(exist(f{i},'file'),sprintf('fail to convert %s (renamed to _tmp%s)',f{i},xt));
    delete(tmp);
end


%%

% func	raw
% afunc	(deletable)
% rfunc	(deletable)
% rafunc	(deletable)
% brafunc	(deletable)
% dbrafunc	min
% kbrafunc	(deletable)
% dkbrafunc	min
% fdbrafunc	(deletable)
% ldbrafunc	(deletable)
% fdkbrafunc	(deletable)
% fldbrafunc	(deletable)
% ldkbrafunc	(deletable)
% v2dbrafunc	alff
% v3dbrafunc	alff
% fldkbrafunc	(deletable)
% v2dkbrafunc	alff
% v2fdbrafunc	conn
% v2ldbrafunc	alff
% v3dkbrafunc	alff
% v3fdbrafunc	conn
% v3ldbrafunc	alff
% wv2dbrafunc	alff
% wv3dbrafunc	alff
% s4v2dbrafunc	alff
% s4v3dbrafunc	alff
% s8v2dbrafunc	alff
% s8v3dbrafunc	alff
% v2fdkbrafunc	conn
% v2fldbrafunc	conn
% v2fldkbrafunc	conn
% v2ldkbrafunc	alff
% v3fdkbrafunc	conn
% v3fldbrafunc	conn
% v3ldkbrafunc	alff
% wv2dkbrafunc	alff
% wv2fdbrafunc	conn
% wv2ldbrafunc	alff
% wv3dkbrafunc	alff
% wv3fdbrafunc	conn
% wv3ldbrafunc	alff
% s4v2dkbrafunc	alff
% s4v2fdbrafunc	conn
% s4v2ldbrafunc	alff
% s4v3dkbrafunc	alff
% s4v3fdbrafunc	conn
% s4v3ldbrafunc	alff
% s8v2dkbrafunc	alff
% s8v2fdbrafunc	conn
% s8v2ldbrafunc	alff
% s8v3dkbrafunc	alff
% s8v3fdbrafunc	conn
% s8v3ldbrafunc	alff
% v3fldkbrafunc	conn
% ws4v2dbrafunc	alff
% ws4v3dbrafunc	alff
% ws8v2dbrafunc	alff
% ws8v3dbrafunc	alff
% wv2fdkbrafunc	conn
% wv2fldbrafunc	conn
% wv2ldkbrafunc	alff
% wv3fdkbrafunc	conn
% wv3fldbrafunc	conn
% wv3ldkbrafunc	alff
% s4v2fdkbrafunc	conn
% s4v2fldbrafunc	conn
% s4v2ldkbrafunc	alff
% s4v3fdkbrafunc	conn
% s4v3fldbrafunc	conn
% s4v3ldkbrafunc	alff
% s8v2fdkbrafunc	conn
% s8v2fldbrafunc	conn
% s8v2ldkbrafunc	alff
% s8v3fdkbrafunc	conn
% s8v3fldbrafunc	conn
% s8v3ldkbrafunc	alff
% ws4v2dkbrafunc	alff
% ws4v2fdbrafunc	conn
% ws4v2ldbrafunc	alff
% ws4v3dkbrafunc	alff
% ws4v3fdbrafunc	conn
% ws4v3ldbrafunc	alff
% ws8v2dkbrafunc	alff
% ws8v2fdbrafunc	conn
% ws8v2ldbrafunc	alff
% ws8v3dkbrafunc	alff
% ws8v3fdbrafunc	conn
% ws8v3ldbrafunc	alff
% wv2fldkbrafunc	conn
% wv3fldkbrafunc	conn
% s4v2fldkbrafunc	conn
% s4v3fldkbrafunc	conn
% s8v2fldkbrafunc	conn
% s8v3fldkbrafunc	conn
% ws4v2fdkbrafunc	conn
% ws4v2fldbrafunc	conn
% ws4v2ldkbrafunc	alff
% ws4v3fdkbrafunc	conn
% ws4v3fldbrafunc	conn
% ws4v3ldkbrafunc	alff
% ws8v2fdkbrafunc	conn
% ws8v2fldbrafunc	conn
% ws8v2ldkbrafunc	alff
% ws8v3fdkbrafunc	conn
% ws8v3fldbrafunc	conn
% ws8v3ldkbrafunc	alff
% ws4v2fldkbrafunc	conn
% ws4v3fldkbrafunc	conn
% ws8v2fldkbrafunc	conn
% ws8v3fldkbrafunc	conn


%%

function xxDiGraph
%%
T = {
    'A'    1      'motion estimation'                        [ 0                   ] '';
    'B'    2      'motion expansion spike'                   [ 1                   ] '';
    'C'    3      'slice timing'                             [ 0                   ] '';
    'D'    4      'realignment'                              [ 3                   ] '';
    'E'    5      'coregistration segmentation'              [ 3                   ] '';
    'G'    6      'deboning'                                 [ 4 5                 ] '';
    'U'    7      'reslicing tissue maps'                    [ 5 6                 ] '';
    'V'    8      'tissue depth'                             [ 7                   ] '';
    'F'    9      'tissue eroding'                           [ 8                   ] '';
    'H'    10     'dvars1 spike'                             [ 6                   ] '';
    'I'    11     'despiking voxelwise'                      [ 6                   ] 'despike';
    'J'    12     'detrending'                               [ 11                  ] 'detrend';
    'J'    13     'detrending'                               [ 6                   ] 'detrend';
    'K'    14     'nuisance extraction'                      [ 9 12                ] 'nuixtr2';
    'K'    15     'nuisance extraction'                      [ 9 13                ] 'nuixtr2';
    'K'    16     'nuisance extraction'                      [ 6 9 12              ] 'nuixtr3';
    'K'    17     'nuisance extraction'                      [ 6 9 13              ] 'nuixtr3';
    'L'    18     'spike interpolation'                      [ 2 6 10 12           ] 'interpspike';
    'L'    19     'spike interpolation'                      [ 2 6 10 13           ] 'interpspike';
    'M'    20     'filtering images'                         [ 18                  ] 'filter';
    'M'    21     'filtering images'                         [ 19                  ] 'filter';
    'M'    22     'filtering images'                         [ 12                  ] 'filter';
    'M'    23     'filtering images'                         [ 13                  ] 'filter';
    'N'    24     'filtering covariates'                     [ 2 14                ] '';
    'N'    25     'filtering covariates'                     [ 2 15                ] '';
    'N'    26     'filtering covariates'                     [ 2 16                ] '';
    'N'    27     'filtering covariates'                     [ 2 17                ] '';
    'O'    28     'regressing out'                           [ 20 24               ] 'regout';
    'O'    29     'regressing out'                           [ 21 25               ] 'regout';
    'O'    30     'regressing out'                           [ 22 24               ] 'regout';
    'O'    31     'regressing out'                           [ 23 25               ] 'regout';
    'O'    32     'regressing out'                           [ 2 18 14             ] 'regout';
    'O'    33     'regressing out'                           [ 2 19 15             ] 'regout';
    'O'    34     'regressing out'                           [ 2 12 14             ] 'regout';
    'O'    35     'regressing out'                           [ 2 13 15             ] 'regout';
    'O'    36     'regressing out'                           [ 20 26               ] 'regout';
    'O'    37     'regressing out'                           [ 21 27               ] 'regout';
    'O'    38     'regressing out'                           [ 22 26               ] 'regout';
    'O'    39     'regressing out'                           [ 23 27               ] 'regout';
    'O'    40     'regressing out'                           [ 2 18 16             ] 'regout';
    'O'    41     'regressing out'                           [ 2 19 17             ] 'regout';
    'O'    42     'regressing out'                           [ 2 12 16             ] 'regout';
    'O'    43     'regressing out'                           [ 2 13 17             ] 'regout';
    'P'    44     'combining spikes to censors'              [ 2 10                ] ''; % force combine censors
    'Q'    45     'smoothing'                                [ 28                  ] 'smooth4';
    'Q'    46     'smoothing'                                [ 29                  ] 'smooth4';
    'Q'    47     'smoothing'                                [ 30                  ] 'smooth4';
    'Q'    48     'smoothing'                                [ 31                  ] 'smooth4';
    'Q'    49     'smoothing'                                [ 32                  ] 'smooth4';
    'Q'    50     'smoothing'                                [ 33                  ] 'smooth4';
    'Q'    51     'smoothing'                                [ 34                  ] 'smooth4';
    'Q'    52     'smoothing'                                [ 35                  ] 'smooth4';
    'Q'    53     'smoothing'                                [ 36                  ] 'smooth4';
    'Q'    54     'smoothing'                                [ 37                  ] 'smooth4';
    'Q'    55     'smoothing'                                [ 38                  ] 'smooth4';
    'Q'    56     'smoothing'                                [ 39                  ] 'smooth4';
    'Q'    57     'smoothing'                                [ 40                  ] 'smooth4';
    'Q'    58     'smoothing'                                [ 41                  ] 'smooth4';
    'Q'    59     'smoothing'                                [ 42                  ] 'smooth4';
    'Q'    60     'smoothing'                                [ 43                  ] 'smooth4';
    'Q'    61     'smoothing'                                [ 28                  ] 'smooth8';
    'Q'    62     'smoothing'                                [ 29                  ] 'smooth8';
    'Q'    63     'smoothing'                                [ 30                  ] 'smooth8';
    'Q'    64     'smoothing'                                [ 31                  ] 'smooth8';
    'Q'    65     'smoothing'                                [ 32                  ] 'smooth8';
    'Q'    66     'smoothing'                                [ 33                  ] 'smooth8';
    'Q'    67     'smoothing'                                [ 34                  ] 'smooth8';
    'Q'    68     'smoothing'                                [ 35                  ] 'smooth8';
    'Q'    69     'smoothing'                                [ 36                  ] 'smooth8';
    'Q'    70     'smoothing'                                [ 37                  ] 'smooth8';
    'Q'    71     'smoothing'                                [ 38                  ] 'smooth8';
    'Q'    72     'smoothing'                                [ 39                  ] 'smooth8';
    'Q'    73     'smoothing'                                [ 40                  ] 'smooth8';
    'Q'    74     'smoothing'                                [ 41                  ] 'smooth8';
    'Q'    75     'smoothing'                                [ 42                  ] 'smooth8';
    'Q'    76     'smoothing'                                [ 43                  ] 'smooth8';
    'R'    77     'normalization'                            [ 5 45                ] 'normalize';
    'R'    78     'normalization'                            [ 5 46                ] 'normalize';
    'R'    79     'normalization'                            [ 5 47                ] 'normalize';
    'R'    80     'normalization'                            [ 5 48                ] 'normalize';
    'R'    81     'normalization'                            [ 5 49                ] 'normalize';
    'R'    82     'normalization'                            [ 5 50                ] 'normalize';
    'R'    83     'normalization'                            [ 5 51                ] 'normalize';
    'R'    84     'normalization'                            [ 5 52                ] 'normalize';
    'R'    85     'normalization'                            [ 5 53                ] 'normalize';
    'R'    86     'normalization'                            [ 5 54                ] 'normalize';
    'R'    87     'normalization'                            [ 5 55                ] 'normalize';
    'R'    88     'normalization'                            [ 5 56                ] 'normalize';
    'R'    89     'normalization'                            [ 5 57                ] 'normalize';
    'R'    90     'normalization'                            [ 5 58                ] 'normalize';
    'R'    91     'normalization'                            [ 5 59                ] 'normalize';
    'R'    92     'normalization'                            [ 5 60                ] 'normalize';
    'R'    93     'normalization'                            [ 5 61                ] 'normalize';
    'R'    94     'normalization'                            [ 5 62                ] 'normalize';
    'R'    95     'normalization'                            [ 5 63                ] 'normalize';
    'R'    96     'normalization'                            [ 5 64                ] 'normalize';
    'R'    97     'normalization'                            [ 5 65                ] 'normalize';
    'R'    98     'normalization'                            [ 5 66                ] 'normalize';
    'R'    99     'normalization'                            [ 5 67                ] 'normalize';
    'R'    100    'normalization'                            [ 5 68                ] 'normalize';
    'R'    101    'normalization'                            [ 5 69                ] 'normalize';
    'R'    102    'normalization'                            [ 5 70                ] 'normalize';
    'R'    103    'normalization'                            [ 5 71                ] 'normalize';
    'R'    104    'normalization'                            [ 5 72                ] 'normalize';
    'R'    105    'normalization'                            [ 5 73                ] 'normalize';
    'R'    106    'normalization'                            [ 5 74                ] 'normalize';
    'R'    107    'normalization'                            [ 5 75                ] 'normalize';
    'R'    108    'normalization'                            [ 5 76                ] 'normalize';
    'R'    109    'normalization'                            [ 5 28                ] 'normalize';
    'R'    110    'normalization'                            [ 5 29                ] 'normalize';
    'R'    111    'normalization'                            [ 5 30                ] 'normalize';
    'R'    112    'normalization'                            [ 5 31                ] 'normalize';
    'R'    113    'normalization'                            [ 5 32                ] 'normalize';
    'R'    114    'normalization'                            [ 5 33                ] 'normalize';
    'R'    115    'normalization'                            [ 5 34                ] 'normalize';
    'R'    116    'normalization'                            [ 5 35                ] 'normalize';
    'R'    117    'normalization'                            [ 5 36                ] 'normalize';
    'R'    118    'normalization'                            [ 5 37                ] 'normalize';
    'R'    119    'normalization'                            [ 5 38                ] 'normalize';
    'R'    120    'normalization'                            [ 5 39                ] 'normalize';
    'R'    121    'normalization'                            [ 5 40                ] 'normalize';
    'R'    122    'normalization'                            [ 5 41                ] 'normalize';
    'R'    123    'normalization'                            [ 5 42                ] 'normalize';
    'R'    124    'normalization'                            [ 5 43                ] 'normalize';
    'R'    125    'normalization'                            [ 5                   ] ; % force normalize structural
    'S'    126    'unwarping atlas'                          [ 5                   ] ;
    'T'    127    'adjacency matrix'                         [ 28 44 126           ] 'adjmat';
    'T'    128    'adjacency matrix'                         [ 29 44 126           ] 'adjmat';
    'T'    129    'adjacency matrix'                         [ 30 44 126           ] 'adjmat';
    'T'    130    'adjacency matrix'                         [ 31 44 126           ] 'adjmat';
    'T'    131    'adjacency matrix'                         [ 32 44 126           ] 'adjmat';
    'T'    132    'adjacency matrix'                         [ 33 44 126           ] 'adjmat';
    'T'    133    'adjacency matrix'                         [ 34 44 126           ] 'adjmat';
    'T'    134    'adjacency matrix'                         [ 35 44 126           ] 'adjmat';
    'T'    135    'adjacency matrix'                         [ 36 44 126           ] 'adjmat';
    'T'    136    'adjacency matrix'                         [ 37 44 126           ] 'adjmat';
    'T'    137    'adjacency matrix'                         [ 38 44 126           ] 'adjmat';
    'T'    138    'adjacency matrix'                         [ 39 44 126           ] 'adjmat';
    'T'    139    'adjacency matrix'                         [ 40 44 126           ] 'adjmat';
    'T'    140    'adjacency matrix'                         [ 41 44 126           ] 'adjmat';
    'T'    141    'adjacency matrix'                         [ 42 44 126           ] 'adjmat';
    'T'    142    'adjacency matrix'                         [ 43 44 126           ] 'adjmat';
    'W'    143    'dvars2'                                   [ 6 12 28             ] '';
    'W'    144    'dvars2'                                   [ 6 13 29             ] '';
    'W'    145    'dvars2'                                   [ 6 12 30             ] '';
    'W'    146    'dvars2'                                   [ 6 13 31             ] '';
    'W'    147    'dvars2'                                   [ 6 12 32             ] '';
    'W'    148    'dvars2'                                   [ 6 13 33             ] '';
    'W'    149    'dvars2'                                   [ 6 12 34             ] '';
    'W'    150    'dvars2'                                   [ 6 13 35             ] '';
    'W'    151    'dvars2'                                   [ 6 12 36             ] '';
    'W'    152    'dvars2'                                   [ 6 13 37             ] '';
    'W'    153    'dvars2'                                   [ 6 12 38             ] '';
    'W'    154    'dvars2'                                   [ 6 13 39             ] '';
    'W'    155    'dvars2'                                   [ 6 12 40             ] '';
    'W'    156    'dvars2'                                   [ 6 13 41             ] '';
    'W'    157    'dvars2'                                   [ 6 12 42             ] '';
    'W'    158    'dvars2'                                   [ 6 13 43             ] '';
    'X'    159    'summary'                                  [ 2 8 10 12 28 44 143 ] 'summary';
    'X'    160    'summary'                                  [ 2 8 10 13 29 44 144 ] 'summary';
    'X'    161    'summary'                                  [ 2 8 10 12 30 44 145 ] 'summary';
    'X'    162    'summary'                                  [ 2 8 10 13 31 44 146 ] 'summary';
    'X'    163    'summary'                                  [ 2 8 10 12 32 44 147 ] 'summary';
    'X'    164    'summary'                                  [ 2 8 10 13 33 44 148 ] 'summary';
    'X'    165    'summary'                                  [ 2 8 10 12 34 44 149 ] 'summary';
    'X'    166    'summary'                                  [ 2 8 10 13 35 44 150 ] 'summary';
    'X'    167    'summary'                                  [ 2 8 10 12 36 44 151 ] 'summary';
    'X'    168    'summary'                                  [ 2 8 10 13 37 44 152 ] 'summary';
    'X'    169    'summary'                                  [ 2 8 10 12 38 44 153 ] 'summary';
    'X'    170    'summary'                                  [ 2 8 10 13 39 44 154 ] 'summary';
    'X'    171    'summary'                                  [ 2 8 10 12 40 44 155 ] 'summary';
    'X'    172    'summary'                                  [ 2 8 10 13 41 44 156 ] 'summary';
    'X'    173    'summary'                                  [ 2 8 10 12 42 44 157 ] 'summary';
    'X'    174    'summary'                                  [ 2 8 10 13 43 44 158 ] 'summary';
    'Y'    175    'tar (atlas_etc anat_etc other_images)'    [ 0                   ] 'tarOtherImg';
    'Y'    176    'tar (Covariates Spikes)'                  [ 0                   ] 'tarNonImg';
    'Y'    177    'tar ((individual funcs))'                 [ 0                   ] 'tarFuncImg';
    }; %#ok<*NBRAK2>
[Label,Source,Detail,Target] = shiDeal(T);
Source = cell2mat(Source);
assert(isequal(Source',1:length(Source)));
Label = shiStrConcat(Label,cellstr(num2str(Source,'%-d')));
Detail = shiStrConcat(Label,':',shiSpace,Detail);
A = false(size(T,1));
for i = 1:length(Source)
    A(Target{i}(Target{i}~=0),Source(i)) = true;
end
G = digraph(A,Label);
plot(G,'NodeLabel',G.Nodes.Name,'Layout','layered');
