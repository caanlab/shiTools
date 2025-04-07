function shiSpmPipe_RestPreproc_WrapperGen

xTr = spm_input('TR (sec)', 1, 'e', [], 1);
while ~isscalar(xTr) || xTr<=0
    xTr = spm_input('[must be number >0] TR (sec)', '+0', 'e', [], 1);
end


%% <Step 1> Slice Timing

spm_input('<Step 1> Slice Timing', 1, 'd');

doSETUP = true;

while doSETUP

    do01_SliceTime = true;
    spm_input('do Slice Timing?', 2, 'b', {'Y'});

    % xSlice_Ta
    % specify the TA parameter in slice timing; usually leave empty for default TR-TR/nSlice for slice indices; ignored for slice times
    MSG = 'TA (sec)';
    xSlice_Ta = spm_input(MSG, [], '[empty]|enter number', [0,1], 1);
    if xSlice_Ta == 0, xSlice_Ta =[]; else
        xSlice_Ta = spm_input(MSG, '+0', 'e', [], 1);
    end

    % xSlice_Order
    % specify slice indices or slice times (ms); leave empty to skip slice-timing and just copy-paste images
    MSG = 'slice order <# or ms>';
    xSlice_Order = spm_input(MSG, [], 'e', [], Inf);

    % xSlice_Ref
    % leave empty to use mid-TR as reference, or specify one slice in the same format as in slice order
    MSG = 'reference slice';
    xSlice_Ref = spm_input(MSG, [], 'mid-TR|enter number', [0,1], 1);
    if xSlice_Ref == 0, xSlice_Ref =[]; else
        xSlice_Ref = spm_input(MSG, '+0', 'e', [], 1);
        while ~isscalar(xSlice_Ref) || xSlice_Ref<=0
            xSlice_Ref = spm_input(['[must be number >=0] ',MSG], '+0', 'e', [], 1);
        end
    end

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end


%% <Step 2> Realignment

spm_input('<Step 2> Realignment', 1, 'd');

doSETUP = true;

while doSETUP

    do02_Realign = spm_input('do Realignment?', 2, 'Y|N', [true,false], 1);
    if ~do02_Realign, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xMotion_AbsMotOption
    % leave empty for default 'Abs_Jenkinson', or specify 'Abs_Power' (see shiSpmPreprocMotCalc)
    MSG = 'absolute motion';
    xMotion_AbsMotOption = spm_input(MSG, [], 'Jenkinson|Power|Other', [1,2,3], 1);
    switch xMotion_AbsMotOption
        case 1, xMotion_AbsMotOption = 'Abs_Jenkinson';
        case 2, xMotion_AbsMotOption = 'Abs_Power';
        case 3, xMotion_AbsMotOption = spm_input(MSG, '+0', 's', '');
        otherwise, xMotion_AbsMotOption = [];
    end

    % xMotion_FdOption
    % leave empty for default 'Fd_Jenkinson', or specify 'Fd_Power' (see shiSpmPreprocMotCalc)
    MSG = 'framewise displacement';
    xMotion_FdOption = spm_input(MSG, [], 'Jenkinson|Power|Other', [1,2,3], 1);
    switch xMotion_FdOption
        case 1, xMotion_FdOption = 'Fd_Jenkinson';
        case 2, xMotion_FdOption = 'Fd_Power';
        case 3, xMotion_FdOption = spm_input(MSG, '+0', 's', '');
        otherwise, xMotion_FdOption = [];
    end

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 3> Skull Stripping

spm_input('<Step 3> Skull Stripping', 1, 'd');

if ~do02_Realign
    doSETUP = false;
    do03_SkullStrip = false;
else
    doSETUP = true;
end

while doSETUP

    do03_SkullStrip = spm_input('do Skull Stripping?', 2, 'Y|N', [true,false], 1);
    if ~do03_SkullStrip, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xDebone_Expr
    % leave empty for default 'i1+i2+i3>0.5' (FSL), or specify 'i1+i2>0.2' (SPM), or specify custom expression (see shiSpmPreprocSkullStrip)
    MSG = 'non-skull mask definition';
    xDebone_Expr = spm_input(MSG, [], 'm', 'GM + WM + CSF > 0.5 (FSL)|GM + WM > 0.2 (SPM)|enter expression', [1,2,3], 1);
    switch xDebone_Expr
        case 1, xDebone_Expr = 'i1+i2+i3>0.5';
        case 2, xDebone_Expr = 'i1+i2>0.2';
        case 3, xDebone_Expr = strrep(strrep(strrep(lower(spm_input(MSG, '+0', 'e', [], 1)),'gm','i1'),'wm','i2'),'csf','i3');
        otherwise, xDebone_Expr = [];
    end
    spm_input(MSG, [], 'b', {xDebone_Expr});

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 4> Despiking

spm_input('<Step 4> Despiking', 1, 'd');

if ~do03_SkullStrip
    doSETUP = false;
    do04_Despike = false;
else
    doSETUP = true;
end

while doSETUP

    do04_Despike = spm_input('do Despiking?', 2, 'Y|N', [true,false], 1);
    if ~do04_Despike, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xDespike_Method 
    % leave empty for default AFNI '3dDespike', or specify 'ArtRepair'
    % xDespike_Parameter 
    % enter {c1,c2,cOrder} for 3dDespike (leave empty for default = {2.5, 4.0, nImg/30}), or {WinSize,Thres} for ArtRepair (leave empty for default = {17, 4}, WinSize must be odd integer)
    MSG = 'despiking method';
    xDespike_Method = spm_input(MSG, [], '3dDespike|ArtRepair', [1,2], 1);
    switch xDespike_Method
        case 1, xDespike_Method = '3dDespike';
            xDespike_Parameter{1} = spm_input('parameter c1:', [], '2.5|enter number', [2.5,NaN], 1);
            if isnan(xDespike_Parameter{1}), xDespike_Parameter{1} = spm_input('parameter c1:', '+0', 'e', 2.5); end
            xDespike_Parameter{2} = spm_input('parameter c2:', [], '4.0|enter number', [4.0,NaN], 1);
            if isnan(xDespike_Parameter{2}), xDespike_Parameter{2} = spm_input('parameter c2:', '+0', 'e', 4.0); end
            xDespike_Parameter{3} = spm_input('parameter cOrder:', [], '#Image/30|enter number', [1,2], 1);
            if xDespike_Parameter{3} == 1, xDespike_Parameter{3} = [];
            elseif xDespike_Parameter{3} == 2, xDespike_Parameter{3} = spm_input('parameter cOrder:', '+0', 'e'); end
        case 2, xDespike_Method = 'ArtRepair';
            xDespike_Parameter{1} = spm_input('parameter WinSize:', [], '17|enter number', [17,NaN], 1);
            if isnan(xDespike_Parameter{1}), xDespike_Parameter{1} = spm_input('parameter WinSize:', '+0', 'e', 17); end
            xDespike_Parameter{2} = spm_input('parameter Threshold:', [], '4|enter number', [4,NaN], 1);
            if isnan(xDespike_Parameter{2}), xDespike_Parameter{2} = spm_input('parameter Threshold:', '+0', 'e', 4); end
        otherwise, xDespike_Method = [];
    end

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 5> Detrending

spm_input('<Step 5> Detrending', 1, 'd');

if ~do03_SkullStrip
    doSETUP = false;
    do05_Detrend = false;
else
    doSETUP = true;
end

while doSETUP

    do05_Detrend = spm_input('do Detrending?', 2, 'Y|N', [true,false], 1);
    if ~do05_Detrend, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xDetrend_Order
    % leave empty for default round(1+xTr*length(func)/150)
    MSG = 'polynomial order';
    xDetrend_Order = spm_input(MSG, [], 'm', '1 + TR * #Image / 150|enter number', [1,2], 1);
    switch xDetrend_Order
        case 1, xDetrend_Order = [];
        case 2, xDetrend_Order = spm_input(MSG, '+0', 'n', [], 1);
        otherwise, xDetrend_Order = [];
    end

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 6> Interpolation

spm_input('<Step 6> Interpolation', 1, 'd');

if ~do05_Detrend
    doSETUP = false;
    do06_Interpolate = false;
else
    doSETUP = true;
end

while doSETUP

    do06_Interpolate = spm_input('do Interpolation?', 2, 'Y|N', [true,false], 1);
    if ~do06_Interpolate, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xDvars_do
    % leave empty for default true
    MSG = 'calculate DVARS?';
    xDvars_do = spm_input(MSG, [], 'Y|N', [true,false], 1);

    if xDvars_do
        % xDvars_SpikeThres
        % leave empty to use FWE-pDVARS<0.05 & delta%D-var>5% (per Afyouni & Nichols, 2018 NeuroImage), or specify an absolute DVARS threshold (e.g. 2) (see shiSpmPreprocDvarsCalc)
        MSG = 'DVARS spike threshold';
        xDvars_SpikeThres = spm_input(MSG, [], 'm', 'pFWE-DVARS < .05 & ∆%D-var > 5%|enter number', [1,2], 1);
        switch xDvars_SpikeThres
            case 1, xDvars_SpikeThres = [];
            case 2, xDvars_SpikeThres = spm_input(MSG, '+0', 'e', [], 1);
            otherwise, xDvars_SpikeThres = [];
        end
    end

    % xMotion_FdSpikeThres
    % leave empty for default 0.5;
    MSG = 'FD spike threshold (mm)';
    xMotion_FdSpikeThres = spm_input(MSG, [], 'br1', '0.5', 0.5);

    % xInterpolate_Method
    % leave empty for default 'Lomb' (Lomb-Scargle periodogram), or specify method for interp1.m ('linear', 'cubic', 'spline', etc)
    MSG = 'Interpolation method';
    xInterpolate_Method = spm_input(MSG, [], 'Lomb|linear|cubic|spline|other', [1,2,3,4,5], 1);
    switch xInterpolate_Method
        case 1, xInterpolate_Method = 'Lomb';
        case 2, xInterpolate_Method = 'linear';
        case 3, xInterpolate_Method = 'cubic';
        case 4, xInterpolate_Method = 'spline';
        case 5, xInterpolate_Method = spm_input(MSG, '+0', 's');
        otherwise, xInterpolate_Method = [];
    end

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 7> Filtering

spm_input('<Step 7> Filtering', 1, 'd');

if ~do05_Detrend
    doSETUP = false;
    do07_Filter = false;
else
    doSETUP = true;
end

while doSETUP

    do07_Filter = spm_input('do Filtering?', 2, 'Y|N', [true,false], 1);
    if ~do07_Filter, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xFilter_HighCutoff
    % leave empty for default 0.01
    MSG = 'high cutoff (Hz)';
    xFilter_HighCutoff = spm_input(MSG, [], 'br1', '0.01', 0.01);

    % xFilter_LowCutoff
    % leave empty for default 0.1
    MSG = 'low cutoff (Hz)';
    xFilter_LowCutoff = spm_input(MSG, [], 'br1', '0.1', 0.1);

    % xFilter_CustCovFiltInd
    % leave empty for default filtering all custom covariates, or 1 * nCov vector of TRUE/FALSE indicating which ones to be filtered
    MSG = 'which custom covariates to skip filter';
    xFilter_CustCovFiltInd = spm_input(MSG, [], 'none|enter indices', [1,2], 1);
    switch xFilter_CustCovFiltInd
        case 1, xFilter_CustCovFiltInd = [];
        case 2, xFilter_CustCovFiltInd = spm_input(MSG, '+0', 'e', [], Inf);
        otherwise, xFilter_CustCovFiltInd = [];
    end

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 8> Regressing Out

spm_input('<Step 8> Regressing Out', 1, 'd');

if ~do05_Detrend
    doSETUP = false;
    do08_RegressOut = false;
else
    doSETUP = true;
end

while doSETUP

    do08_RegressOut = spm_input('do Regressing Out?', 2, 'Y|N', [true,false], 1);
    if ~do08_RegressOut, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xErode_Keep
    % leave empty for default 0.1, i.e. keep 10% of volume, or specify custom cutoff (see shiSpmPreprocErode)
    MSG = 'deep tissue threshold';
    xErode_Keep = spm_input(MSG, [], '10% deepest|enter number', [1,2], 1);
    switch xErode_Keep
        case 1, xErode_Keep = 0.1;
        case 2, xErode_Keep = spm_input([MSG,' eg .1 for top 10%'], '+0', 'r', [], 1);
        otherwise, xErode_Keep = [];
    end

    % nui2 or nui3?
    MSG = 'nuisance signals';
    yNui = spm_input(MSG, [], 'm', 'WM+CSF|WB+WM+CSF|both WM+CSF & WB+WM+CSF', [1,2,3], 1);

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 9> Smoothing

spm_input('<Step 9> Smoothing', 1, 'd');

if ~do08_RegressOut
    doSETUP = false;
    do09_Smooth = false;
else
    doSETUP = true;
end

while doSETUP

    do09_Smooth = spm_input('do Smoothing?', 2, 'Y|N', [true,false], 1);
    if ~do09_Smooth, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % fwhm4 or fwhm8?
    MSG = 'full width at half maximum (FWHM)';
    yFwhm = spm_input(MSG, [], 'm', '4mm|8mm|both 4mm & 8mm', [1,2,3], 3);

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Step 10> Normalization

spm_input('<Step 10> Normalization', 1, 'd');

if ~do08_RegressOut
    doSETUP = false;
    do10_Normalise = false;
else
    doSETUP = true;
end

while doSETUP

    do10_Normalise = spm_input('do Normalization?', 2, 'Y|N', [true,false], 1);
    if ~do10_Normalise, doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]); continue; end

    % xNormalize_VoxelSize
    % leave empty for default [3 3 3] (mm)
    MSG = 'normalized voxel size (mm)';
    xNormalize_VoxelSize = spm_input(MSG, [], 'br1', '1|2|3', [1,2,3], 3);

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% <Final Step> Adjacency Matrix

spm_input('<Final Step> Adjacency Matrix', 1, 'd');

if ~do08_RegressOut
    doSETUP = false;
else
    doSETUP = true;
end

while doSETUP

    % xAdjMat_XtrSummFunc
    % leave empty for default 'mean', or specify 'med' or 'eig'
    MSG = 'parcel summary method';
    xAdjMat_XtrSummFunc = spm_input(MSG, [], 'mean|median|eigen', [1,2,3], 1);
    switch xAdjMat_XtrSummFunc
        case 1, xAdjMat_XtrSummFunc = 'mean';
        case 2, xAdjMat_XtrSummFunc = 'med';
        case 3, xAdjMat_XtrSummFunc = 'eig';
        otherwise, xAdjMat_XtrSummFunc = [];
    end

    % xAdjMat_CorrMethod
    % leave empty for default 'pearson' (same as @(x)atanh(corr(x,'type','pearson')), or specify similarly 'spearman', or any function handle
    MSG = 'connectivity method';
    xAdjMat_CorrMethod = spm_input(MSG, [], 'm', 'Pearson -> Fisher-Z|Spearman -> Fisher-Z|enter function handle', [1,2,3], 1);
    switch xAdjMat_CorrMethod
        case 1, xAdjMat_CorrMethod = 'pearson';
        case 2, xAdjMat_CorrMethod = 'spearman';
        case 3, xAdjMat_CorrMethod = spm_input(MSG, '+0', 'e', '@(x)atanh(corr(x))');
        otherwise, xAdjMat_CorrMethod = [];
    end

    atlas_ALL = 'AAL|AAL_90|AAL3|AAL3_140|HCPMMP|HOA|Power264|Neuromorphometrics_Brain|Schaefer1000Parcel07Net|Schaefer1000Parcel17Net|Yeo7NetLiberal|Yeo17NetLiberal|all';
    atlas = spm_input('choose atlas', [], 'm', atlas_ALL, 1:13, 13);
    atlas_ALL = split(atlas_ALL,'|');
    atlas = atlas_ALL(atlas);
    if isequal(atlas,'all')
        atlas = atlas_ALL; %#ok<NASGU>
    end

    doSETUP = spm_input(' ', '+2', 'Next|Reset', [false,true]);

end

%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

spm_input('Generating wrapper script...', 1, 'd');
pause(0.1);

FIN = cellstr([
    shiIf(do10_Normalise   , 'w', '');
    shiIf(do09_Smooth      , 's', '');
    shiIf(do08_RegressOut  , 'v', '');
    shiIf(do07_Filter      , 'f', '');
    shiIf(do06_Interpolate , 'l', '');
    shiIf(do05_Detrend     , 'd', '');
    shiIf(do04_Despike     , 'k', '');
    shiIf(do03_SkullStrip  , 'b', '');
    shiIf(do02_Realign     , 'r', '');
    shiIf(do01_SliceTime   , 'a', '');
    ]');
if exist('yNui','var')
    switch yNui
        case 1, FIN = cellstr(shiStrRepl(FIN, 'v', 'v2', false));
        case 2, FIN = cellstr(shiStrRepl(FIN, 'v', 'v3', false));
        otherwise, FIN = [shiStrRepl(FIN, 'v', 'v2', false); shiStrRepl(FIN, 'v', 'v3', false)];
    end
end
if exist('yFwhm','var')
    switch yFwhm
        case 1, FIN = cellstr(shiStrRepl(FIN, 's', 's4', false));
        case 2, FIN = cellstr(shiStrRepl(FIN, 's', 's8', false));
        otherwise, FIN = [cellstr(shiStrRepl(FIN, 's', 's4', false)); cellstr(shiStrRepl(FIN, 's', 's8', false))];
    end
end
if ~exist('atlas','var')
    atlas = []; %#ok<NASGU>
end

[xStepInd,xStepInd_letter,xStepInd_detail] = deal(cell(size(FIN)));
for i = 1:length(FIN)
    [xStepInd{i},xStepInd_letter{i},xStepInd_detail{i}] = shipipe_get_steps(FIN{i});
end
xStepInd = cat(1,xStepInd{:});
xStepInd_letter = cat(1,xStepInd_letter{:});
xStepInd_detail = cat(1,xStepInd_detail{:});
[xStepInd, idx] = unique(xStepInd);
xStepInd_letter = xStepInd_letter(idx);
xStepInd_detail = xStepInd_detail(idx);

clear doSETUP MSG idx atlas_ALL

%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

[WRAP, WRAP_PATH] = uiputfile({'*.m'},'Save wrapper script','pipe_wrapper.m');
if isequal(WRAP,0), error('aborted'); end
WRAP = fullfile(WRAP_PATH,WRAP);

try xPARAM.xTr                    = xTr                    ; catch, warning('unspecified parameter: %s', 'xTr                   '); xPARAM.xTr                    = []; end
try xPARAM.xSlice_Ta              = xSlice_Ta              ; catch, warning('unspecified parameter: %s', 'xSlice_Ta             '); xPARAM.xSlice_Ta              = []; end
try xPARAM.xSlice_Order           = xSlice_Order           ; catch, warning('unspecified parameter: %s', 'xSlice_Order          '); xPARAM.xSlice_Order           = []; end
try xPARAM.xSlice_Ref             = xSlice_Ref             ; catch, warning('unspecified parameter: %s', 'xSlice_Ref            '); xPARAM.xSlice_Ref             = []; end
try xPARAM.xDvars_do              = xDvars_do              ; catch, warning('unspecified parameter: %s', 'xDvars_do             '); xPARAM.xDvars_do              = []; end
try xPARAM.xDespike_Method        = xDespike_Method        ; catch, warning('unspecified parameter: %s', 'xDespike_Method       '); xPARAM.xDespike_Method        = []; end
try xPARAM.xInterpolate_Method    = xInterpolate_Method    ; catch, warning('unspecified parameter: %s', 'xInterpolate_Method   '); xPARAM.xInterpolate_Method    = []; end
try xPARAM.xNormalize_VoxelSize   = xNormalize_VoxelSize   ; catch, warning('unspecified parameter: %s', 'xNormalize_VoxelSize  '); xPARAM.xNormalize_VoxelSize   = []; end
try xPARAM.xMotion_AbsMotOption   = xMotion_AbsMotOption   ; catch, warning('unspecified parameter: %s', 'xMotion_AbsMotOption  '); xPARAM.xMotion_AbsMotOption   = []; end
try xPARAM.xMotion_FdOption       = xMotion_FdOption       ; catch, warning('unspecified parameter: %s', 'xMotion_FdOption      '); xPARAM.xMotion_FdOption       = []; end
try xPARAM.xMotion_FdSpikeThres   = xMotion_FdSpikeThres   ; catch, warning('unspecified parameter: %s', 'xMotion_FdSpikeThres  '); xPARAM.xMotion_FdSpikeThres   = []; end
try xPARAM.xDebone_Expr           = xDebone_Expr           ; catch, warning('unspecified parameter: %s', 'xDebone_Expr          '); xPARAM.xDebone_Expr           = []; end
try xPARAM.xErode_Keep            = xErode_Keep            ; catch, warning('unspecified parameter: %s', 'xErode_Keep           '); xPARAM.xErode_Keep            = []; end
try xPARAM.xDvars_SpikeThres      = xDvars_SpikeThres      ; catch, warning('unspecified parameter: %s', 'xDvars_SpikeThres     '); xPARAM.xDvars_SpikeThres      = []; end
try xPARAM.xDespike_Parameter     = xDespike_Parameter     ; catch, warning('unspecified parameter: %s', 'xDespike_Parameter    '); xPARAM.xDespike_Parameter     = []; end
try xPARAM.xDetrend_Order         = xDetrend_Order         ; catch, warning('unspecified parameter: %s', 'xDetrend_Order        '); xPARAM.xDetrend_Order         = []; end
try xPARAM.xFilter_HighCutoff     = xFilter_HighCutoff     ; catch, warning('unspecified parameter: %s', 'xFilter_HighCutoff    '); xPARAM.xFilter_HighCutoff     = []; end
try xPARAM.xFilter_LowCutoff      = xFilter_LowCutoff      ; catch, warning('unspecified parameter: %s', 'xFilter_LowCutoff     '); xPARAM.xFilter_LowCutoff      = []; end
try xPARAM.xFilter_CustCovFiltInd = xFilter_CustCovFiltInd ; catch, warning('unspecified parameter: %s', 'xFilter_CustCovFiltInd'); xPARAM.xFilter_CustCovFiltInd = []; end
try xPARAM.xAdjMat_XtrSummFunc    = xAdjMat_XtrSummFunc    ; catch, warning('unspecified parameter: %s', 'xAdjMat_XtrSummFunc   '); xPARAM.xAdjMat_XtrSummFunc    = []; end
try xPARAM.xAdjMat_CorrMethod     = xAdjMat_CorrMethod     ; catch, warning('unspecified parameter: %s', 'xAdjMat_CorrMethod    '); xPARAM.xAdjMat_CorrMethod     = []; end %#ok<STRNU>

matlab.io.saveVariablesToScript(WRAP, {'xStepInd','atlas','xPARAM'});

WRAP_VAR = shiTxtRead(WRAP);

WRAP_TXT0 = {
    'anat = ;'
    ''
    'func = ;'
    ''
    'CustomCov = []; % optional'
    'CustomSpike = []; % optional'
    ''
    };

WRAP_TXT2 = [
    {''}
    shiStrConcat('%', shiSpace, cellstr(num2str(xStepInd,'[%d]')), shiSpace, '[', xStepInd_letter, ']', shiSpace, xStepInd_detail)
    {''}
    ];

WRAP_TXT3 = {
    'shiSpmPipe_RestPreproc(xStepInd,anat,func,atlas,CustomCov,CustomSpike,xPARAM);'
    };

WRAP_ALL = [WRAP_TXT0; WRAP_VAR; WRAP_TXT2; WRAP_TXT3];

fid = fopen(WRAP,'w');
fprintf(fid, '%s\n', WRAP_ALL{:});
fclose(fid);

spm_input('Done.', 2, 'd');


%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function [FIN_req,FIN_req_letter,FIN_req_detail] = shipipe_get_steps(FIN)

% e.g., FIN = 'ws4v2fldkbra';

[NUM, NUM2, DETAIL, IN,OUT] = shiDeal({
    1,    'A',  'motion estimation',                                      {'FUNC'},                                                                                                                                                                                                                                                  {'Rp_func1'};
    2,    'B',  'calculating motion parameters and motion-spikes',        {'Rp_func1'},                                                                                                                                                                                                                                              {'Mot24_Rp_func1', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'SpikeFd_Rp_func1'};
    3,    'C',  'slice timing',                                           {'FUNC'},                                                                                                                                                                                                                                                  {'afunc'};
    4,    'D',  'realignment',                                            {'afunc'},                                                                                                                                                                                                                                                 {'rafunc', 'Mean_afunc1'};
    5,    'E',  'coregistration and segmentation',                        {'ANAT', 'Mean_afunc1'},                                                                                                                                                                                                                                   {'y_anat', 'iy_anat', 'c1anat', 'c2anat', 'c3anat', 'manat'};
    6,    'G',  'skull stripping',                                        {'rafunc', 'c1anat', 'c2anat', 'c3anat'},                                                                                                                                                                                                                  {'brafunc', 'MaskDebone_rafunc1'};
    7,    'U',  'reslicing tissue maps',                                  {'brafunc', 'c1anat', 'c2anat', 'c3anat'},                                                                                                                                                                                                                 {'Resliced_c1anat', 'Resliced_c2anat', 'Resliced_c3anat'};
    8,    'V',  'calculating tissue depth',                               {'Resliced_c1anat', 'Resliced_c2anat', 'Resliced_c3anat'},                                                                                                                                                                                                 {'TisDep_Resliced_c1anat'};
    9,    'F',  'tissue eroding',                                         {'TisDep_Resliced_c1anat'},                                                                                                                                                                                                                                {'ec2_TisDep_Resliced_c1anat', 'ec3_TisDep_Resliced_c1anat'};
    10,   'H',  'calculating initial DVARS and DVARS-spikes',             {'brafunc', 'MaskDebone_rafunc1'},                                                                                                                                                                                                                         {'Dvars_brafunc1', 'SpikeDvars_brafunc1'};
    11,   'I',  'voxelwise despiking',                                    {'brafunc', 'MaskDebone_rafunc1'},                                                                                                                                                                                                                         {'kbrafunc'};
    12,   'J',  'detrending',                                             {'kbrafunc', 'FUNC'},                                                                                                                                                                                                                                      {'dkbrafunc', 'Mean_kbrafunc1'};
    13,   'J',  'detrending',                                             {'brafunc', 'FUNC'},                                                                                                                                                                                                                                       {'dbrafunc', 'Mean_brafunc1'};
    14,   'K',  'nuisance extraction',                                    {'dkbrafunc', 'ec2_TisDep_Resliced_c1anat', 'ec3_TisDep_Resliced_c1anat', 'Nui2_dkbrafunc1'},                                                                                                                                                              {'Nui2_dkbrafunc1'};
    15,   'K',  'nuisance extraction',                                    {'dbrafunc', 'ec2_TisDep_Resliced_c1anat', 'ec3_TisDep_Resliced_c1anat', 'Nui2_dbrafunc1'},                                                                                                                                                                {'Nui2_dbrafunc1'};
    16,   'K',  'nuisance extraction',                                    {'dkbrafunc', 'ec2_TisDep_Resliced_c1anat', 'ec3_TisDep_Resliced_c1anat', 'MaskDebone_rafunc1', 'Nui3_dkbrafunc1'},                                                                                                                                        {'Nui3_dkbrafunc1'};
    17,   'K',  'nuisance extraction',                                    {'dbrafunc', 'ec2_TisDep_Resliced_c1anat', 'ec3_TisDep_Resliced_c1anat', 'MaskDebone_rafunc1', 'Nui3_dbrafunc1'},                                                                                                                                          {'Nui3_dbrafunc1'};
    18,   'L',  'interpolating spike volumes',                            {'dkbrafunc', 'MaskDebone_rafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE'},                                                                                                                                                             {'ldkbrafunc'};
    19,   'L',  'interpolating spike volumes',                            {'dbrafunc', 'MaskDebone_rafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE'},                                                                                                                                                              {'ldbrafunc'};
    20,   'M',  'filtering images',                                       {'ldkbrafunc'},                                                                                                                                                                                                                                            {'fldkbrafunc'};
    21,   'M',  'filtering images',                                       {'ldbrafunc'},                                                                                                                                                                                                                                             {'fldbrafunc'};
    22,   'M',  'filtering images',                                       {'dkbrafunc'},                                                                                                                                                                                                                                             {'fdkbrafunc'};
    23,   'M',  'filtering images',                                       {'dbrafunc'},                                                                                                                                                                                                                                              {'fdbrafunc'};
    24,   'N',  'filtering nuisance brain signals',                       {'Nui2_dkbrafunc1'},                                                                                                                                                                                                                                       {'fNui2_dkbrafunc1'};
    25,   'N',  'filtering nuisance brain signals',                       {'Nui2_dbrafunc1'},                                                                                                                                                                                                                                        {'fNui2_dbrafunc1'};
    26,   'N',  'filtering nuisance brain signals',                       {'Nui3_dkbrafunc1'},                                                                                                                                                                                                                                       {'fNui3_dkbrafunc1'};
    27,   'N',  'filtering nuisance brain signals',                       {'Nui3_dbrafunc1'},                                                                                                                                                                                                                                        {'fNui3_dbrafunc1'};
    28,   'N',  'filtering motion parameters and other covariates',       {'Mot24_Rp_func1'},                                                                                                                                                                                                                                        {'fMot24_Rp_func1'};
    29,   'N',  'filtering custom covariates',                            {'CUSTOMCOV'},                                                                                                                                                                                                                                             {'fCUSTOMCOV'};
    30,   'O',  'regressing out covariates',                              {'fldkbrafunc', 'fMot24_Rp_func1', 'fNui2_dkbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                      {'v2fldkbrafunc', 'v2Spm_fldkbrafunc1'};
    31,   'O',  'regressing out covariates',                              {'fldbrafunc', 'fMot24_Rp_func1', 'fNui2_dbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                        {'v2fldbrafunc', 'v2Spm_fldbrafunc1'};
    32,   'O',  'regressing out covariates',                              {'fdkbrafunc', 'fMot24_Rp_func1', 'fNui2_dkbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                       {'v2fdkbrafunc', 'v2Spm_fdkbrafunc1'};
    33,   'O',  'regressing out covariates',                              {'fdbrafunc', 'fMot24_Rp_func1', 'fNui2_dbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                         {'v2fdbrafunc', 'v2Spm_fdbrafunc1'};
    34,   'O',  'regressing out covariates',                              {'ldkbrafunc', 'Mot24_Rp_func1', 'Nui2_dkbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                          {'v2ldkbrafunc', 'v2Spm_ldkbrafunc1'};
    35,   'O',  'regressing out covariates',                              {'ldbrafunc', 'Mot24_Rp_func1', 'Nui2_dbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                            {'v2ldbrafunc', 'v2Spm_ldbrafunc1'};
    36,   'O',  'regressing out covariates',                              {'dkbrafunc', 'Mot24_Rp_func1', 'Nui2_dkbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                           {'v2dkbrafunc', 'v2Spm_dkbrafunc1'};
    37,   'O',  'regressing out covariates',                              {'dbrafunc', 'Mot24_Rp_func1', 'Nui2_dbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                             {'v2dbrafunc', 'v2Spm_dbrafunc1'};
    38,   'O',  'regressing out covariates',                              {'fldkbrafunc', 'fMot24_Rp_func1', 'fNui3_dkbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                      {'v3fldkbrafunc', 'v3Spm_fldkbrafunc1'};
    39,   'O',  'regressing out covariates',                              {'fldbrafunc', 'fMot24_Rp_func1', 'fNui3_dbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                        {'v3fldbrafunc', 'v3Spm_fldbrafunc1'};
    40,   'O',  'regressing out covariates',                              {'fdkbrafunc', 'fMot24_Rp_func1', 'fNui3_dkbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                       {'v3fdkbrafunc', 'v3Spm_fdkbrafunc1'};
    41,   'O',  'regressing out covariates',                              {'fdbrafunc', 'fMot24_Rp_func1', 'fNui3_dbrafunc1', 'fCUSTOMCOV'},                                                                                                                                                                                         {'v3fdbrafunc', 'v3Spm_fdbrafunc1'};
    42,   'O',  'regressing out covariates',                              {'ldkbrafunc', 'Mot24_Rp_func1', 'Nui3_dkbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                          {'v3ldkbrafunc', 'v3Spm_ldkbrafunc1'};
    43,   'O',  'regressing out covariates',                              {'ldbrafunc', 'Mot24_Rp_func1', 'Nui3_dbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                            {'v3ldbrafunc', 'v3Spm_ldbrafunc1'};
    44,   'O',  'regressing out covariates',                              {'dkbrafunc', 'Mot24_Rp_func1', 'Nui3_dkbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                           {'v3dkbrafunc', 'v3Spm_dkbrafunc1'};
    45,   'O',  'regressing out covariates',                              {'dbrafunc', 'Mot24_Rp_func1', 'Nui3_dbrafunc1', 'CUSTOMCOV'},                                                                                                                                                                                             {'v3dbrafunc', 'v3Spm_dbrafunc1'};
    46,   'P',  'combining spikes to censors',                            {'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE'},                                                                                                                                                                                                {'Censor_brafunc1'};
    47,   'Q',  'smoothing',                                              {'v2fldkbrafunc'},                                                                                                                                                                                                                                         {'s4v2fldkbrafunc'};
    48,   'Q',  'smoothing',                                              {'v2fldbrafunc'},                                                                                                                                                                                                                                          {'s4v2fldbrafunc'};
    49,   'Q',  'smoothing',                                              {'v2fdkbrafunc'},                                                                                                                                                                                                                                          {'s4v2fdkbrafunc'};
    50,   'Q',  'smoothing',                                              {'v2fdbrafunc'},                                                                                                                                                                                                                                           {'s4v2fdbrafunc'};
    51,   'Q',  'smoothing',                                              {'v2ldkbrafunc'},                                                                                                                                                                                                                                          {'s4v2ldkbrafunc'};
    52,   'Q',  'smoothing',                                              {'v2ldbrafunc'},                                                                                                                                                                                                                                           {'s4v2ldbrafunc'};
    53,   'Q',  'smoothing',                                              {'v2dkbrafunc'},                                                                                                                                                                                                                                           {'s4v2dkbrafunc'};
    54,   'Q',  'smoothing',                                              {'v2dbrafunc'},                                                                                                                                                                                                                                            {'s4v2dbrafunc'};
    55,   'Q',  'smoothing',                                              {'v3fldkbrafunc'},                                                                                                                                                                                                                                         {'s4v3fldkbrafunc'};
    56,   'Q',  'smoothing',                                              {'v3fldbrafunc'},                                                                                                                                                                                                                                          {'s4v3fldbrafunc'};
    57,   'Q',  'smoothing',                                              {'v3fdkbrafunc'},                                                                                                                                                                                                                                          {'s4v3fdkbrafunc'};
    58,   'Q',  'smoothing',                                              {'v3fdbrafunc'},                                                                                                                                                                                                                                           {'s4v3fdbrafunc'};
    59,   'Q',  'smoothing',                                              {'v3ldkbrafunc'},                                                                                                                                                                                                                                          {'s4v3ldkbrafunc'};
    60,   'Q',  'smoothing',                                              {'v3ldbrafunc'},                                                                                                                                                                                                                                           {'s4v3ldbrafunc'};
    61,   'Q',  'smoothing',                                              {'v3dkbrafunc'},                                                                                                                                                                                                                                           {'s4v3dkbrafunc'};
    62,   'Q',  'smoothing',                                              {'v3dbrafunc'},                                                                                                                                                                                                                                            {'s4v3dbrafunc'};
    63,   'Q',  'smoothing',                                              {'v2fldkbrafunc'},                                                                                                                                                                                                                                         {'s8v2fldkbrafunc'};
    64,   'Q',  'smoothing',                                              {'v2fldbrafunc'},                                                                                                                                                                                                                                          {'s8v2fldbrafunc'};
    65,   'Q',  'smoothing',                                              {'v2fdkbrafunc'},                                                                                                                                                                                                                                          {'s8v2fdkbrafunc'};
    66,   'Q',  'smoothing',                                              {'v2fdbrafunc'},                                                                                                                                                                                                                                           {'s8v2fdbrafunc'};
    67,   'Q',  'smoothing',                                              {'v2ldkbrafunc'},                                                                                                                                                                                                                                          {'s8v2ldkbrafunc'};
    68,   'Q',  'smoothing',                                              {'v2ldbrafunc'},                                                                                                                                                                                                                                           {'s8v2ldbrafunc'};
    69,   'Q',  'smoothing',                                              {'v2dkbrafunc'},                                                                                                                                                                                                                                           {'s8v2dkbrafunc'};
    70,   'Q',  'smoothing',                                              {'v2dbrafunc'},                                                                                                                                                                                                                                            {'s8v2dbrafunc'};
    71,   'Q',  'smoothing',                                              {'v3fldkbrafunc'},                                                                                                                                                                                                                                         {'s8v3fldkbrafunc'};
    72,   'Q',  'smoothing',                                              {'v3fldbrafunc'},                                                                                                                                                                                                                                          {'s8v3fldbrafunc'};
    73,   'Q',  'smoothing',                                              {'v3fdkbrafunc'},                                                                                                                                                                                                                                          {'s8v3fdkbrafunc'};
    74,   'Q',  'smoothing',                                              {'v3fdbrafunc'},                                                                                                                                                                                                                                           {'s8v3fdbrafunc'};
    75,   'Q',  'smoothing',                                              {'v3ldkbrafunc'},                                                                                                                                                                                                                                          {'s8v3ldkbrafunc'};
    76,   'Q',  'smoothing',                                              {'v3ldbrafunc'},                                                                                                                                                                                                                                           {'s8v3ldbrafunc'};
    77,   'Q',  'smoothing',                                              {'v3dkbrafunc'},                                                                                                                                                                                                                                           {'s8v3dkbrafunc'};
    78,   'Q',  'smoothing',                                              {'v3dbrafunc'},                                                                                                                                                                                                                                            {'s8v3dbrafunc'};
    79,   'R',  'normalization',                                          {'s4v2fldkbrafunc', 'y_anat'},                                                                                                                                                                                                                             {'ws4v2fldkbrafunc'};
    80,   'R',  'normalization',                                          {'s4v2fldbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws4v2fldbrafunc'};
    81,   'R',  'normalization',                                          {'s4v2fdkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws4v2fdkbrafunc'};
    82,   'R',  'normalization',                                          {'s4v2fdbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws4v2fdbrafunc'};
    83,   'R',  'normalization',                                          {'s4v2ldkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws4v2ldkbrafunc'};
    84,   'R',  'normalization',                                          {'s4v2ldbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws4v2ldbrafunc'};
    85,   'R',  'normalization',                                          {'s4v2dkbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws4v2dkbrafunc'};
    86,   'R',  'normalization',                                          {'s4v2dbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'ws4v2dbrafunc'};
    87,   'R',  'normalization',                                          {'s4v3fldkbrafunc', 'y_anat'},                                                                                                                                                                                                                             {'ws4v3fldkbrafunc'};
    88,   'R',  'normalization',                                          {'s4v3fldbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws4v3fldbrafunc'};
    89,   'R',  'normalization',                                          {'s4v3fdkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws4v3fdkbrafunc'};
    90,   'R',  'normalization',                                          {'s4v3fdbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws4v3fdbrafunc'};
    91,   'R',  'normalization',                                          {'s4v3ldkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws4v3ldkbrafunc'};
    92,   'R',  'normalization',                                          {'s4v3ldbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws4v3ldbrafunc'};
    93,   'R',  'normalization',                                          {'s4v3dkbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws4v3dkbrafunc'};
    94,   'R',  'normalization',                                          {'s4v3dbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'ws4v3dbrafunc'};
    95,   'R',  'normalization',                                          {'s8v2fldkbrafunc', 'y_anat'},                                                                                                                                                                                                                             {'ws8v2fldkbrafunc'};
    96,   'R',  'normalization',                                          {'s8v2fldbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws8v2fldbrafunc'};
    97,   'R',  'normalization',                                          {'s8v2fdkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws8v2fdkbrafunc'};
    98,   'R',  'normalization',                                          {'s8v2fdbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws8v2fdbrafunc'};
    99,   'R',  'normalization',                                          {'s8v2ldkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws8v2ldkbrafunc'};
    100,  'R',  'normalization',                                          {'s8v2ldbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws8v2ldbrafunc'};
    101,  'R',  'normalization',                                          {'s8v2dkbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws8v2dkbrafunc'};
    102,  'R',  'normalization',                                          {'s8v2dbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'ws8v2dbrafunc'};
    103,  'R',  'normalization',                                          {'s8v3fldkbrafunc', 'y_anat'},                                                                                                                                                                                                                             {'ws8v3fldkbrafunc'};
    104,  'R',  'normalization',                                          {'s8v3fldbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws8v3fldbrafunc'};
    105,  'R',  'normalization',                                          {'s8v3fdkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws8v3fdkbrafunc'};
    106,  'R',  'normalization',                                          {'s8v3fdbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws8v3fdbrafunc'};
    107,  'R',  'normalization',                                          {'s8v3ldkbrafunc', 'y_anat'},                                                                                                                                                                                                                              {'ws8v3ldkbrafunc'};
    108,  'R',  'normalization',                                          {'s8v3ldbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws8v3ldbrafunc'};
    109,  'R',  'normalization',                                          {'s8v3dkbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'ws8v3dkbrafunc'};
    110,  'R',  'normalization',                                          {'s8v3dbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'ws8v3dbrafunc'};
    111,  'R',  'normalization',                                          {'v2fldkbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'wv2fldkbrafunc'};
    112,  'R',  'normalization',                                          {'v2fldbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'wv2fldbrafunc'};
    113,  'R',  'normalization',                                          {'v2fdkbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'wv2fdkbrafunc'};
    114,  'R',  'normalization',                                          {'v2fdbrafunc', 'y_anat'},                                                                                                                                                                                                                                 {'wv2fdbrafunc'};
    115,  'R',  'normalization',                                          {'v2ldkbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'wv2ldkbrafunc'};
    116,  'R',  'normalization',                                          {'v2ldbrafunc', 'y_anat'},                                                                                                                                                                                                                                 {'wv2ldbrafunc'};
    117,  'R',  'normalization',                                          {'v2dkbrafunc', 'y_anat'},                                                                                                                                                                                                                                 {'wv2dkbrafunc'};
    118,  'R',  'normalization',                                          {'v2dbrafunc', 'y_anat'},                                                                                                                                                                                                                                  {'wv2dbrafunc'};
    119,  'R',  'normalization',                                          {'v3fldkbrafunc', 'y_anat'},                                                                                                                                                                                                                               {'wv3fldkbrafunc'};
    120,  'R',  'normalization',                                          {'v3fldbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'wv3fldbrafunc'};
    121,  'R',  'normalization',                                          {'v3fdkbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'wv3fdkbrafunc'};
    122,  'R',  'normalization',                                          {'v3fdbrafunc', 'y_anat'},                                                                                                                                                                                                                                 {'wv3fdbrafunc'};
    123,  'R',  'normalization',                                          {'v3ldkbrafunc', 'y_anat'},                                                                                                                                                                                                                                {'wv3ldkbrafunc'};
    124,  'R',  'normalization',                                          {'v3ldbrafunc', 'y_anat'},                                                                                                                                                                                                                                 {'wv3ldbrafunc'};
    125,  'R',  'normalization',                                          {'v3dkbrafunc', 'y_anat'},                                                                                                                                                                                                                                 {'wv3dkbrafunc'};
    126,  'R',  'normalization',                                          {'v3dbrafunc', 'y_anat'},                                                                                                                                                                                                                                  {'wv3dbrafunc'};
    127,  'R',  'normalization of structural images',                     {'ANAT', 'manat', 'y_anat'},                                                                                                                                                                                                                               {'wanat', 'wmanat'};
    128,  'S',  'unwarping atlas',                                        {'ATLAS', 'iy_anat'},                                                                                                                                                                                                                                      {'uatlas'};
    129,  'T',  'calculating adjacency matrix',                           {'v2fldkbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v2fldkbrafunc1_uatlas'},                                                                                                                                                                {'Adj_v2fldkbrafunc1_uatlas'};
    130,  'T',  'calculating adjacency matrix',                           {'v2fldbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v2fldbrafunc1_uatlas'},                                                                                                                                                                  {'Adj_v2fldbrafunc1_uatlas'};
    131,  'T',  'calculating adjacency matrix',                           {'v2fdkbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v2fdkbrafunc1_uatlas'},                                                                                                                                                                  {'Adj_v2fdkbrafunc1_uatlas'};
    132,  'T',  'calculating adjacency matrix',                           {'v2fdbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v2fdbrafunc1_uatlas'},                                                                                                                                                                    {'Adj_v2fdbrafunc1_uatlas'};
    133,  'T',  'calculating adjacency matrix',                           {'v3fldkbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v3fldkbrafunc1_uatlas'},                                                                                                                                                                {'Adj_v3fldkbrafunc1_uatlas'};
    134,  'T',  'calculating adjacency matrix',                           {'v3fldbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v3fldbrafunc1_uatlas'},                                                                                                                                                                  {'Adj_v3fldbrafunc1_uatlas'};
    135,  'T',  'calculating adjacency matrix',                           {'v3fdkbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v3fdkbrafunc1_uatlas'},                                                                                                                                                                  {'Adj_v3fdkbrafunc1_uatlas'};
    136,  'T',  'calculating adjacency matrix',                           {'v3fdbrafunc', 'uatlas', 'ATLAS_LABEL', 'Censor_brafunc1', 'Adj_v3fdbrafunc1_uatlas'},                                                                                                                                                                    {'Adj_v3fdbrafunc1_uatlas'};
    137,  'W',  'calculating final DVARS',                                {'v2fldkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                 {'Dvars_v2fldkbrafunc1'};
    138,  'W',  'calculating final DVARS',                                {'v2fldbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                   {'Dvars_v2fldbrafunc1'};
    139,  'W',  'calculating final DVARS',                                {'v2fdkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                  {'Dvars_v2fdkbrafunc1'};
    140,  'W',  'calculating final DVARS',                                {'v2fdbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                    {'Dvars_v2fdbrafunc1'};
    141,  'W',  'calculating final DVARS',                                {'v2ldkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                  {'Dvars_v2ldkbrafunc1'};
    142,  'W',  'calculating final DVARS',                                {'v2ldbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                    {'Dvars_v2ldbrafunc1'};
    143,  'W',  'calculating final DVARS',                                {'v2dkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                   {'Dvars_v2dkbrafunc1'};
    144,  'W',  'calculating final DVARS',                                {'v2dbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                     {'Dvars_v2dbrafunc1'};
    145,  'W',  'calculating final DVARS',                                {'v3fldkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                 {'Dvars_v3fldkbrafunc1'};
    146,  'W',  'calculating final DVARS',                                {'v3fldbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                   {'Dvars_v3fldbrafunc1'};
    147,  'W',  'calculating final DVARS',                                {'v3fdkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                  {'Dvars_v3fdkbrafunc1'};
    148,  'W',  'calculating final DVARS',                                {'v3fdbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                    {'Dvars_v3fdbrafunc1'};
    149,  'W',  'calculating final DVARS',                                {'v3ldkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                  {'Dvars_v3ldkbrafunc1'};
    150,  'W',  'calculating final DVARS',                                {'v3ldbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                    {'Dvars_v3ldbrafunc1'};
    151,  'W',  'calculating final DVARS',                                {'v3dkbrafunc', 'MaskDebone_rafunc1', 'Mean_kbrafunc1'},                                                                                                                                                                                                   {'Dvars_v3dkbrafunc1'};
    152,  'W',  'calculating final DVARS',                                {'v3dbrafunc', 'MaskDebone_rafunc1', 'Mean_brafunc1'},                                                                                                                                                                                                     {'Dvars_v3dbrafunc1'};
    153,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v2fldkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2fldkbrafunc1', 'CUSTOMCOV', 'v2Spm_fldkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},      {'PreprocSumm_v2fldkbrafunc1'};
    154,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v2fldbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2fldbrafunc1', 'CUSTOMCOV', 'v2Spm_fldbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},          {'PreprocSumm_v2fldbrafunc1'};
    155,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v2fdkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2fdkbrafunc1', 'CUSTOMCOV', 'v2Spm_fdkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},         {'PreprocSumm_v2fdkbrafunc1'};
    156,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v2fdbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2fdbrafunc1', 'CUSTOMCOV', 'v2Spm_fdbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},             {'PreprocSumm_v2fdbrafunc1'};
    157,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v2ldkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2ldkbrafunc1', 'CUSTOMCOV', 'v2Spm_ldkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},         {'PreprocSumm_v2ldkbrafunc1'};
    158,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v2ldbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2ldbrafunc1', 'CUSTOMCOV', 'v2Spm_ldbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},             {'PreprocSumm_v2ldbrafunc1'};
    159,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v2dkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2dkbrafunc1', 'CUSTOMCOV', 'v2Spm_dkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},            {'PreprocSumm_v2dkbrafunc1'};
    160,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v2dbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v2dbrafunc1', 'CUSTOMCOV', 'v2Spm_dbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},                {'PreprocSumm_v2dbrafunc1'};
    161,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v3fldkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3fldkbrafunc1', 'CUSTOMCOV', 'v3Spm_fldkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},      {'PreprocSumm_v3fldkbrafunc1'};
    162,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v3fldbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3fldbrafunc1', 'CUSTOMCOV', 'v3Spm_fldbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},          {'PreprocSumm_v3fldbrafunc1'};
    163,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v3fdkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3fdkbrafunc1', 'CUSTOMCOV', 'v3Spm_fdkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},         {'PreprocSumm_v3fdkbrafunc1'};
    164,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v3fdbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3fdbrafunc1', 'CUSTOMCOV', 'v3Spm_fdbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},             {'PreprocSumm_v3fdbrafunc1'};
    165,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v3ldkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3ldkbrafunc1', 'CUSTOMCOV', 'v3Spm_ldkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},         {'PreprocSumm_v3ldkbrafunc1'};
    166,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v3ldbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3ldbrafunc1', 'CUSTOMCOV', 'v3Spm_ldbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},             {'PreprocSumm_v3ldbrafunc1'};
    167,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dkbrafunc', 'v3dkbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3dkbrafunc1', 'CUSTOMCOV', 'v3Spm_dkbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},            {'PreprocSumm_v3dkbrafunc1'};
    168,  'X',  'summary',                                                {'TisDep_Resliced_c1anat', 'dbrafunc', 'v3dbrafunc', 'AbsMot_Rp_func1', 'Fd_Rp_func1', 'Dvars_brafunc1', 'Dvars_v3dbrafunc1', 'CUSTOMCOV', 'v3Spm_dbrafunc1', 'SpikeFd_Rp_func1', 'SpikeDvars_brafunc1', 'CUSTOMSPIKE', 'Censor_brafunc1'},                {'PreprocSumm_v3dbrafunc1'};
    });

assert(isequal(cell2mat(NUM),(1:length(NUM))'));

REQ = cell(size(IN));
for i = 1:length(REQ)
    REQ{i} = i;
    in = IN{i};
    req = unique([REQ{i}, shipipe_get_in_from_out(in,OUT)]);
    while ~isequal(req, REQ{i})
        in = [IN{setdiff(req, REQ{i})}];
        REQ{i} = req;
        req = unique([REQ{i}, shipipe_get_in_from_out(in,OUT)]);
    end
end

YIE = cell(size(OUT));
for i = 1:length(OUT)
    tmp = OUT{i}(endsWith(OUT{i},'func') | startsWith(OUT{i},'Adj_') | startsWith(OUT{i},'Dvars_') | startsWith(OUT{i},'PreprocSumm_'));
    if isempty(tmp)
        YIE(i) = {''};
    else
        YIE(i) = tmp;
    end
end
YIE = shiStrRepl(YIE,'func1','',false);
YIE = shiStrRepl(YIE,'func','',false);

FIN_index = find(matches(YIE,FIN));
assert(isscalar(FIN_index));

if contains(FIN,'w'), FIN_index = [FIN_index; find(strcmp(DETAIL,'normalization of structural images'))]; end % need to do structural normalization

FIN0 = FIN;
FIN0 = shiStrRepl(FIN0,'w','',false);
FIN0 = shiStrRepl(FIN0,'s4','',false);
FIN0 = shiStrRepl(FIN0,'s8','',false);

FIN_index = [FIN_index; find(matches(YIE,['Adj_',FIN0,'_uatlas']))];
FIN_index = [FIN_index; find(matches(YIE,['Dvars_',FIN0]))];
FIN_index = [FIN_index; find(matches(YIE,['PreprocSumm_',FIN0]))];

FIN_req = [REQ{FIN_index'}]';
FIN_req = unique(FIN_req);

FIN_req_letter = NUM2(FIN_req);
FIN_req_detail = DETAIL(FIN_req);

function req = shipipe_get_in_from_out(in,OUT)
req = nan(size(in));
for i = 1:length(in)
    tmp = find(cellfun(@(x)ismember(in(i),x),OUT));
    req(i) = shiIf(isempty(tmp),NaN,tmp);
end
req = unique(req);
req(isnan(req)) = [];

