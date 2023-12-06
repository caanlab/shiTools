function [TabDat,TabDat_PRINT,xSPM] = shiSpmResultTable_spm(PathSpm,ContrastIndex,Mask,isInclusiveMasking,PVal,PCorrection,K,PeakDis,PeakNumber,AtlasName)

% reads SPM results and returns activation table
%
% PathSpm            = string for path of SPM.mat
% ContrastIndex      = index of contrast, e.g. 3 for spmT_0003 or spmF_0003
% Mask               = mask image file(s), or contrast indices, or empty []
% isInclusiveMasking = true for inclusive (default), false for exclusive
% PVal               = voxel p value
% PCorrection        = 'FWE' or 'FDR' or 'none'
% K                  = cluster size
% PeakDis            = minimal distance between peaks (mm) (default=8, SPM_default=8)
% PeakNumber         = maximal number of peaks returned (default=Inf, SPM_default=3)
%
% ##########
% Zhenhao Shi  2023-10-19

if ~exist('Mask','var') || isempty(Mask)
    Mask = [];
elseif ischar(Mask)
    Mask = cellstr(char(Mask));
elseif isnumeric(Mask)
    Mask = double(Mask(:));
end

if ~exist('isInclusiveMasking','var') || isempty(isInclusiveMasking)
    isInclusiveMasking = true;
end

if ~exist('PCorrection','var') || isempty(PCorrection) || ~ischar(PCorrection)
    xCorr = 'none';
else
    switch lower(PCorrection)
        case 'none'
            xCorr = 'none';
        case 'fwe'
            xCorr = 'FWE';
        case 'fdr'
            xCorr = 'FDR';
        otherwise
            xCorr = 'none';
    end
end

if ~exist('PeakDis','var') || isempty(PeakDis) || PeakDis<=0
    PeakDis = 8;
end

if ~exist('PeakNumber','var') || isempty(PeakNumber) || PeakNumber<=0
    PeakNumber = Inf;
end

if ~exist('AtlasName','var')
    AtlasName = '';
end

PWD = pwd;
spm('defaults','FMRI')


try
    xSPM.swd       = PathSpm;
    xSPM.Ic        = ContrastIndex;
    xSPM.Im        = Mask;
    xSPM.Ex        = ~isInclusiveMasking;
    xSPM.u         = PVal;
    xSPM.k         = K;
    xSPM.thresDesc = xCorr;
    xSPM.n         = 1; % conjunction
    [~,xSPM] = spm_getSPM(xSPM);
catch
    [~,xSPM] = spm_results_ui('Setup');
end

TabDat = spm_list('Table',xSPM,PeakNumber,PeakDis);

TabDat = [TabDat.hdr(1:2,3:end);TabDat.dat(:,3:end)];

if exist('AtlasName','var') && ~isempty(AtlasName)
    if isempty(strcmp(AtlasName,shiSpmAnatLabel('avail')))
        AtlasList = shiSpmAnatLabel('avail');
        AtlasName = AtlasList{spm_input('Select an atlas','+1','m',[AtlasList{1},sprintf('|%s',AtlasList{2:end})])};
    end
    Lab = shiSpmAnatLabel(AtlasName,cellfun(@(x)(x(:)'),TabDat(3:end,10),'UniformOutput',false),5,true);
    TabDat = [TabDat,[{'Label';''};Lab]];
end

TabDat_PRINT = [TabDat,shiSpmResultTable_TablePrintText(TabDat)];

fprintf('\n%s\n',fullfile(xSPM.swd,xSPM.Vspm.fname));
if iscellstr(Mask) %#ok<ISCLSTR>
    MaskInfo = sprintf('%s  ',Mask{:});
elseif isnumeric(Mask) && ~isempty(Mask)
    MaskInfo = ['Contrast(s) -',sprintf(' %d',Mask)];
else
    MaskInfo = 'none';
end
fprintf('\nMask           : %s (%s)',MaskInfo,shiIf(isInclusiveMasking,'inclusive','exclusive'));
fprintf('\nThreshold      : %s > %.2f, p(%s) < %g',xSPM.STATstr,xSPM.u,xCorr,PVal);
fprintf('\nExtent         : K >= %d',K);
fprintf('\nConnectedness  : 18');
fprintf('\nPeak distance  : %d mm apart\n\n',PeakDis);
if exist('AtlasName','var') && ~isempty(AtlasName)
    fprintf('Atlas          : %s\n',AtlasName);
end

cd(PWD)


function TABLE_PRINT_TEXT = shiSpmResultTable_TablePrintText(TABLE_PRINT)

if size(TABLE_PRINT,2) ~= 11
    TABLE_PRINT_TEXT = {};
    return;
end

xK = 3;
xP = 9;
xCoord = 10;
xLab = 11;

TABLE_PRINT_TEXT = cell(size(TABLE_PRINT,1),1);
TABLE_PRINT_TEXT{1} = 'Text';
for i = 3:size(TABLE_PRINT,1)
    if ~isempty(TABLE_PRINT{i,xK})
        TABLE_PRINT_TEXT{i,1} = sprintf('%s (k=%d, Z=%.2f, x/y/z=%d/%d/%d)', TABLE_PRINT{i,xLab}, TABLE_PRINT{i,xK}, spm_invNcdf(1-TABLE_PRINT{i,xP}), TABLE_PRINT{i,xCoord}(1), TABLE_PRINT{i,xCoord}(2), TABLE_PRINT{i,xCoord}(3));
    else
        TABLE_PRINT_TEXT{i,1} = sprintf('%s (Z=%.2f, x/y/z=%d/%d/%d)', TABLE_PRINT{i,xLab}, spm_invNcdf(1-TABLE_PRINT{i,xP}), TABLE_PRINT{i,xCoord}(1), TABLE_PRINT{i,xCoord}(2), TABLE_PRINT{i,xCoord}(3));
    end
end

