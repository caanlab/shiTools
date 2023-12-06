function [TABLE_SPM,TABLE_PRINT,Vo,Yo,YL,ParamStruct] = shiSpmResultTable(Image,Mask,Statistic,Df1,Df2,Sign,VoxThres,ClustThres,Conn,PeakDistance,PeakNumber,AtlasName)

% reads nifti statistical map and returns activation results table
%
% Image = statistic map filename, e.g. 'spmT_0001.img'
% Mask = .img/.nii filenames or ''
% Statistic = 'T' or 'F' or 'X'
% Df1 = df for t-test or df_1 for f-test;
% Df2 = [] for t-test or df_2 for f-test
% Sign = 'pos' or 'neg' or 'both' for t-test, ignored for f-test
% VoxThres = voxel threshold (p-value or t/F statistic, or arbitrary cutoff)
% ClustThres = cluster size
% Conn = connectedness for cluster forming - 6, 18 or 26; SPM default is 18
% PeakDistance = minimal distance between peaks (mm); SPM default is 8
% PeakNumber = max number of peaks per cluster; default = Inf
% AtlasName = one of the atlas provided in shiTools, see shiSpmAnatLabel; default = 'Neuromorphometrics'
% 
% TABLE_SPM format:
%       col.  1 - LnNum
%       col.  2 - ClustRk
%       col.  3 - NumPeak
%       col.  4 - Extent
%       col.  5 - VoxRk
%       col.  6 - Stat
%       col.  7 - Z
%       col.  8 - MNIx
%       col.  9 - MNIy
%       col. 10 - MNIz
% Vo    - Image output header
% Yo    - Thresholded image value (to save thresholded map: spm_write_vol(Vo,Yo))
% YL    - N-nary image value (to save n-nary map: spm_write_vol(Vo,YL))
%
% ##########
% Zhenhao Shi  2020-03-12
% ##########

if nargin == 0
    
    
    
    Image = spm_select(1,'image','Select an image file');
    [x1,x2,x3,~] = spm_fileparts(Image);
    Image = fullfile(x1,[x2,x3]);

    isMask = spm_input('Inclusive mask : ',-1, 'b', '50% Brain|40% Grey|Other|None',[1,2,3,4],1);
    switch isMask
        case 1
            Mask = which('shi_c123_50.nii');%,',1'];
        case 2
            Mask = which('shi_c1_40.nii');%,',1'];
        case 3
            Mask = spm_select(1,'image','Select a mask file');
            if isequal(Mask(end-1:end),',1')
                Mask = Mask(1:end-2);
            end
        case 4
            Mask = '';
    end
    
    Stat_spm_input = spm_input('Test (T/F/X) : ','+1', 'b', 'T|F|X',[1,2,3],1);
    Statistic = shiIf(Stat_spm_input==1,'T',shiIf(Stat_spm_input==2,'F','X'));

    Df1 = -1;
    Df2 = -1;
    if strcmpi(Statistic,'T')
        try
            FF=spm_vol(Image);
            FF=textscan(FF.descrip,'SPM{T_[%f');
            Df1 = FF{1};
            if ~(Df1>0)
                error('x');
            end
            Df1 = spm_input('T-test df : ','+1', 'n',Df1);
        catch
            Df1 = spm_input('T-test df : ','+1', 'n');
        end
    elseif strcmpi(Statistic,'F')
        try
            FF=spm_vol(Image);
            FF=textscan(FF.descrip,'SPM{F_[%f,%f');
            Df1 = FF{1};
            Df2 = FF{2};
            if ~(Df1>0) || ~(Df2>0)
                error('x');
            end
            Df1 = spm_input('F-test df_1 : ','+1', 'n',Df1);
            Df2 = spm_input('F-test df_2 : ','+1', 'n',Df2);
        catch
            Df1 = spm_input('F-test df_1 : ','+1', 'n',1);
            Df2 = spm_input('F-test df_2 : ','+1', 'n');
        end
    end

    Sign = '';
    if strcmpi(Statistic,'T') || strcmpi(Statistic,'X')
        Sign_Opt = {'pos','neg','both'};
        Sign = Sign_Opt{spm_input('Sign (pos/neg/both) :','+1', 'b', 'pos|neg|both',[1,2,3],1)};
    end

    VoxThres = spm_input('Threshold : ','+1', 'r', 0.001);

    ClustThres = spm_input('Cluster size : ','+1', 'n', 50);

    Conn = spm_input('Connectedness (6/18/26): ','+1', 'r', 18);

    PeakDistance = spm_input('Peak distance : ','+1', 'n', 8);
    PeakNumber = spm_input('Peak number : ','+1', 'n', Inf);
    
    AtlasList = shiSpmAnatLabel('avail');
    AtlasName = AtlasList{spm_input('Select an atlas','+1','m',[AtlasList{1},sprintf('|%s',AtlasList{2:end})])};
   
    fprintf('\nCOMMAND:\n    [T_SPM,T_PRINT,V,Y] = shiSpmResultTable(''%s'',''%s'',''%s'',%d,%d,''%s'',%g,%d,%d,%d,%d,''%s'');\n\n',Image,Mask,Statistic,Df1,Df2,Sign,VoxThres,ClustThres,Conn,PeakDistance,PeakNumber,AtlasName);
    
    
    
else



    if ~exist('Image','var') || isempty(Image) || ~exist(Image,'file')
        Image = spm_select(1,'image','Select an image file');
        [x1,x2,x3,~] = spm_fileparts(Image);
        Image = fullfile(x1,[x2,x3]);
    end
    if ~spm_existfile(Image)
        error('must choose an image');
    end

    if ~exist('Mask','var') || ( ~exist(Mask,'file') && ~isempty(Mask) )
        isMask = spm_input('Inclusive mask : ',-1, 'b', '50% Brain|40% Grey|Other|None',[1,2,3,4],1);
        switch isMask
            case 1
                Mask = which('shi_c123_50.nii');%,',1'];
            case 2
                Mask = which('shi_c1_40.nii');%,',1'];
            case 3
                Mask = spm_select(1,'image','Select a mask file');
                if isequal(a.mask(end-1:end),',1')
                    Mask = a.mask(1:end-2);
                end
            case 4
                Mask = '';
        end
    end

    if ~exist('Statistic','var') || isempty(Statistic)
        Statistic=spm_vol(Image);
        Statistic=textscan(Statistic.descrip,'SPM{%c');
        Statistic=Statistic{1};
        if ~isequal(Statistic,'T') && ~isequal(Statistic,'F')
            Statistic='X';
        end
    end

    if strcmpi(Statistic,'T') && (~exist('Df1','var') || isempty(Df1))
        Df2=-1;
        try
            FF=spm_vol(Image);
            FF=textscan(FF.descrip,'SPM{T_[%f');
            Df1 = FF{1};
            if isempty(Df1) || ~(Df1>0)
                error('x');
            end
        catch
            Df1 = spm_input('T-test df : ','+1', 'n');
        end
    elseif strcmpi(Statistic,'F') && (~exist('Df2','var') || isempty(Df1) || isempty(Df2))
        try
            FF=spm_vol(Image);
            FF=textscan(FF.descrip,'SPM{F_[%f,%f');
            Df1 = FF{1};
            Df2 = FF{2};
            if isempty(Df1) || isempty(Df2) ||~(Df1>0) || ~(Df2>0)
                error('x');
            end
        catch
            Df1 = spm_input('F-test df_1 : ','+1', 'n',1);
            Df2 = spm_input('F-test df_2 : ','+1', 'n');
        end
    elseif strcmpi(Statistic,'X')
        Df1 = -1;
        Df2 = -1;
    end

    if ( strcmpi(Statistic,'T') || strcmpi(Statistic,'X') ) && (~exist('Sign','var') || (~strcmpi(Sign(1:3),'pos') && ~strcmpi(Sign(1:3),'neg') && ~strcmpi(Sign(1:3),'bot')))
        Sign = 'both';
    end

    if ~exist('VoxThres','var') || isempty(VoxThres)
        VoxThres = spm_input('Threshold : ','+1', 'r', 0.001);
    end

    if ~exist('ClustThres','var') || isempty(ClustThres) ||  ~(ClustThres >= 0)
        ClustThres = spm_input('Cluster size : ','+1', 'n', 50);
    end

    if ~exist('Conn','var') || isempty(Conn) || (Conn-6)*(Conn-18)*(Conn-26) ~= 0
        Conn = 18;
    end

    if ~exist('PeakDistance','var') || isempty(PeakDistance) || ~(PeakDistance >= 0)
        PeakDistance = 8;
    end
    
    if ~exist('PeakNumber','var') || isempty(PeakNumber) || ~(PeakNumber >= 0)
        PeakNumber = Inf;
    end
    
    if ~exist('AtlasName','var')
        AtlasName = 'Neuromorphometrics';
    end

    
    
end





fprintf('\n%s\n',char(shiFullFileName(Image)));


%% mask

V = spm_vol(Image);

if isempty(Mask)
    [Y,XYZ] = spm_read_vols(V);
    fprintf('\nNo mask\n');
elseif ~exist(Mask,'file')
    fprintf(2,'\nCANNOT FIND MASK FILE %s\n',Mask);
    [Y,XYZ] = spm_read_vols(V);
else
    [Y,XYZ] = shiSpmResultTable_mask(Mask, Image);
    Y = Y{1};
    fprintf('\nMask         : %s\n',Mask);
end


%% critical threshold value

if ( strcmpi(Statistic,'T') || strcmpi(Statistic,'X') ) && strcmpi(Sign(1:3),'bot')
    [TABLE_SPM1,TABLE_PRINT1,Vo,Yo1,YL1,ParamStruct] = shiSpmResultTable(Image,Mask,Statistic,Df1,Df2,'pos',VoxThres,ClustThres,Conn,PeakDistance,PeakNumber,AtlasName);
    [TABLE_SPM2,TABLE_PRINT2,~  ,Yo2,YL2] = shiSpmResultTable(Image,Mask,Statistic,Df1,Df2,'neg',VoxThres,ClustThres,Conn,PeakDistance,PeakNumber,AtlasName);
    TABLE_SPM = [TABLE_SPM1;TABLE_SPM2];
    TABLE_PRINT = [TABLE_PRINT1;TABLE_PRINT2];
    Yo1(isnan(Yo1)) = 0;
    Yo2(isnan(Yo2)) = 0;
    Yo = nan(size(Yo1));
    Yo(YL1>0) = Yo1(YL1>0);
    Yo(YL2>0) = Yo2(YL2>0);
    YL = YL1-YL2;
    ParamStruct.Sign = 'both';
    return;
end

if ( strcmpi(Statistic,'T') || strcmpi(Statistic,'X') ) && strcmpi(Sign(1:3),'neg')
    Y = -Y;
end

if strcmpi(Statistic,'T')
    if VoxThres>1
        Threshold = VoxThres;
        VoxThres = 1-spm_Tcdf(VoxThres,Df1);
    else
        Threshold = spm_invTcdf(1-VoxThres,Df1);
    end
    fprintf('Threshold      : %st(%d) > %.2f, p < %g\n',shiIf(strcmpi(Sign(1:3),'neg'),'-','+'),Df1,Threshold,VoxThres);
elseif strcmpi(Statistic,'F')
    if VoxThres>1
        Threshold = VoxThres;
        VoxThres = 1-spm_Fcdf(VoxThres,Df1,Df2);
    else
        Threshold = spm_invFcdf(1-VoxThres,Df1,Df2);
    end
    fprintf('Threshold      : F(%d,%d) > %.2f, p < %g\n',Df1,Df2,Threshold,VoxThres);
elseif strcmpi(Statistic,'X')
%     if strcmpi(Sign(1:3),'neg')
%         Threshold = -Pval;
%     else
        Threshold = VoxThres;
%     end
    fprintf('Threshold      : %sX > %.2f\n',shiIf(strcmpi(Sign(1:3),'neg'),'-','+'),Threshold);
end

fprintf('Extent         : K >= %d\n',ClustThres);
fprintf('Connectedness  : %d\n',Conn);
fprintf('Peaks          : %d mm apart\n',PeakDistance);


%% binary and cluster image

BW = (Y > Threshold) .* 1;
[YL,nCluster] = spm_bwlabel(BW,Conn);

[YL,nCluster,ClusterSize] = shiSpmResultTable_ClusterSize(YL,nCluster,ClustThres);

Yo = nan(size(Y));
Yo(YL>0) = Y(YL>0);
Vo = V;
[pp,nn,xx] = fileparts(V.fname);
Vo.fname = fullfile(pp,[nn,'_thres',xx]);
Vo = rmfield(Vo,'private');


fprintf('# clusters     : %d\n',nCluster);


%% peak coordinates and values

xExtent = cell(nCluster,1);
xStat = cell(nCluster,1);
xXYZ = cell(nCluster,1);
xZ = cell(nCluster,1);
xNumPeak = cell(nCluster,1);
xClustRk = cell(nCluster,1);
xVoxRk = cell(nCluster,1);

xThePeak = zeros(nCluster,1);

for i = 1:nCluster

    Y_cluster = Y .* (YL==i);
    Y_cluster(YL~=i) = NaN;
    idx = spm_get_lm(Y_cluster);
%     numel(idx)

    xXYZ{i} = XYZ(:,idx);
    xStat{i} = Y_cluster(idx);
    [xStat{i},IX] = sort(xStat{i},'descend');
    xXYZ{i} = xXYZ{i}(:,IX);

%     length(xStatistic{i})

    [xStat{i},xXYZ{i}] = shiSpmResultTable_PeakApart(xStat{i},xXYZ{i},PeakDistance);
    nPeak = length(xStat{i});
    nPeak = min(nPeak,PeakNumber);
    xStat{i} = xStat{i}(1:nPeak)';
    xXYZ{i} = xXYZ{i}(:,1:nPeak)';

    xExtent{i} = repmat(ClusterSize(i),nPeak,1);
    xZ{i} = shiSpmResultTable_Stat2Z(xStat{i},Statistic,Df1,Df2);
    xNumPeak{i} = repmat(nPeak,nPeak,1);
    xVoxRk{i} = transpose(1:nPeak);
    xClustRk{i} = repmat(i,nPeak,1);
    
    xThePeak(i) = xStat{i}(1);

end

[~,IX] = sort(xThePeak,'descend');
xXYZ = xXYZ(IX);
xStat = xStat(IX);
xExtent = xExtent(IX);
xZ = xZ(IX);
xNumPeak = xNumPeak(IX);
xVoxRk = xVoxRk(IX);
xClustRk = xClustRk(IX);

fprintf('# peaks        : %d\n',numel(cell2mat(xNumPeak)));


%% table

TABLE_NUMonly = cell2mat([xClustRk,xNumPeak,xExtent,xVoxRk,xStat,xZ,xXYZ]);
TABLE_NUMonly = [transpose(1:size(TABLE_NUMonly,1)),TABLE_NUMonly];
TABLE_PRINT = [{'LnNum','ClustRk','NumPeak','Extent','VoxRk','Stat','Z','MNIx','MNIy','MNIz'};num2cell(TABLE_NUMonly)];
TABLE_SPM = TABLE_PRINT;
for i = size(TABLE_SPM,1):-1:3
    if TABLE_SPM{i,2}==TABLE_SPM{i-1,2}
        TABLE_SPM(i,2:4) = {[],[],[]};
    end
end
TABLE_SPM = TABLE_SPM(:,[1,4,6,7,8,9,10]);

if exist('AtlasName','var') && ~isempty(AtlasName)
    if isempty(strcmp(AtlasName,shiSpmAnatLabel('avail')))
        AtlasList = shiSpmAnatLabel('avail');
        AtlasName = AtlasList{spm_input('Select an atlas','+1','m',[AtlasList{1},sprintf('|%s',AtlasList{2:end})])};
    end
    fprintf('Atlas          : %s\n',AtlasName);
    Lab = shiSpmAnatLabel(AtlasName,cell2mat(xXYZ),5,true);
    TABLE_PRINT = [TABLE_PRINT,[{'Label'};Lab]];
    TABLE_PRINT = [TABLE_PRINT,shiSpmResultTable_TablePrintText(TABLE_PRINT)];
    TABLE_SPM = [TABLE_SPM,[{'Label'};Lab]];
else
    AtlasName = '';
end

if strcmpi(Sign(1:3),'neg')
    Yo = -Yo;
end


%% struct

ParamStruct.Image           = Image;
ParamStruct.Mask            = Mask;
ParamStruct.Statistic       = Statistic;
ParamStruct.Df1             = Df1;
ParamStruct.Df2             = Df2;
ParamStruct.Sign            = Sign;
ParamStruct.VoxThres        = Threshold;
ParamStruct.ClustThres      = ClustThres;
ParamStruct.Conn            = Conn;
ParamStruct.PeakDistance    = PeakDistance;
ParamStruct.PeakNumber      = PeakNumber;
ParamStruct.AtlasName       = AtlasName;



function [L2,nCluster2,ClusterSize2] = shiSpmResultTable_ClusterSize(L,nCluster,K)

L2 = zeros(size(L));
ClusterSize = nan(nCluster,1);

for i = 1:nCluster
    ClusterSize(i) = numel(L(L==i));
end

[ClusterSize2,IX] = sort(ClusterSize,'descend');
IX = IX(ClusterSize2>=K);
ClusterSize2 = ClusterSize2(ClusterSize2>=K);
nCluster2 = length(ClusterSize2);


for i = 1:nCluster2
    L2 = L2 + (L==IX(i)) .* i;
end


function [PeakY,PeakXYZ] = shiSpmResultTable_PeakApart(PeakY,PeakXYZ,PeakDistance)

i = 2;

while i <= length(PeakY)

    XYZ_prev = PeakXYZ(:,1:(i-1));
    XYZ_curr = repmat(PeakXYZ(:,i),1,size(XYZ_prev,2));

    Dist = sqrt(sum((XYZ_prev - XYZ_curr) .^ 2));

    if any(Dist <= PeakDistance)
        PeakY = PeakY([1:i-1,i+1:end]);
        PeakXYZ = PeakXYZ(:,[1:i-1,i+1:end]);
    else
        i = i + 1;
    end

end


function Z = shiSpmResultTable_Stat2Z(Stat,Type,Df1,Df2)
switch lower(Type)
    case 't'
        Z = spm_t2z(Stat,Df1);
    case 'f'
        Z = spm_invNcdf(spm_Fcdf(Stat,Df1,Df2));
    case 'x'
        Z = nan(size(Stat));
end        


function [Y,XYZ] = shiSpmResultTable_mask(Masks, Images)
% Mask images, into space that is same as Images{1} or Images(1,:)
% spm_mask.m 4419 2011-08-03 18:42:35Z guillaume $

%-Parameters & Arguments
%--------------------------------------------------------------------------

if iscell(Masks)
    Masks = char(Masks(:));
end
if iscell(Images)
    Images = char(Images(:));
end

V1 = spm_vol(Masks);
V2 = spm_vol(Images);
[~,XYZ] = spm_read_vols(V2(1));

thresh = zeros(numel(V1),1);

m1 = numel(V1);
m2 = numel(V2);

%-Create headers
%--------------------------------------------------------------------------
VO = V2;
Y = cell(m2,1);
for i=1:m2
    [pth,nm,ext,num] = spm_fileparts(deblank(VO(i).fname));
    VO(i).fname      = fullfile(pth,['m', nm, ext, num]);
    VO(i).descrip    = 'Masked';
    VO(i).mat        = VO(1).mat;
    VO(i).dim(1:3)   = VO(1).dim(1:3);
    Y{i} = nan(VO(1).dim(1:3));
end
% VO  = spm_create_vol(VO);
M   = VO(1).mat;
dim = VO(1).dim(1:3);

%-Compute masked images
%--------------------------------------------------------------------------

for j=1:dim(3)

    msk = true(dim(1:2));
    Mi  = spm_matrix([0 0 j]);

    % Load slice j from all images
    for i=1:m1
        M1  = M\V1(i).mat\Mi;
        %if sum((M1(:)-Mi(:)).^2<eps) M1 = Mi; end;

        img = spm_slice_vol(V1(i),M1,dim(1:2),[0 NaN]);

        msk = msk & isfinite(img);

        if ~spm_type(V1(i).dt(1),'nanrep')
            msk = msk & (img ~= 0);
        end

        msk = msk & (img > thresh(i));
    end

    % Write the images
    for i=1:m2
        M1        = M\V2(i).mat\Mi;
        img       = spm_slice_vol(V2(i),M1,dim(1:2),[1 0]);
        img(~msk) = NaN;
        Y{i}(:,:,j) = img;
%         VO(i)     = spm_write_plane(VO(i),img,j);
    end
end


function TABLE_PRINT_TEXT = shiSpmResultTable_TablePrintText(TABLE_PRINT)

if size(TABLE_PRINT,2) ~= 11
    TABLE_PRINT_TEXT = {};
    return;
end

TABLE_PRINT_TEXT = cell(size(TABLE_PRINT,1),1);
for i = 1:size(TABLE_PRINT,1)
    if ischar(TABLE_PRINT{i,1})
        TABLE_PRINT_TEXT{i,1} = 'Text';
    elseif ~isequal(TABLE_PRINT{i,2},TABLE_PRINT{i-1,2})
        TABLE_PRINT_TEXT{i,1} = sprintf('%s (k=%d, Z=%.2f, x/y/z=%d/%d/%d)', TABLE_PRINT{i,11}, TABLE_PRINT{i,4}, TABLE_PRINT{i,7}, TABLE_PRINT{i,8}, TABLE_PRINT{i,9}, TABLE_PRINT{i,10});
    else
        TABLE_PRINT_TEXT{i,1} = sprintf('%s (Z=%.2f, x/y/z=%d/%d/%d)', TABLE_PRINT{i,11}, TABLE_PRINT{i,7}, TABLE_PRINT{i,8}, TABLE_PRINT{i,9}, TABLE_PRINT{i,10});
    end
end



