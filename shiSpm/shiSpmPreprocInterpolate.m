function outImg = shiSpmPreprocInterpolate(Img,Mask,Flag,Tr,InterpMethod,Prefix,existAction)
%
% interpolates flagged timepoints
%
%  Img          : cellstr image names
%  Mask         : char mask image (default = [])
%  Flag         : logical matrix, or .txt file(s) of such, that flag timepoints to be interpolated (TRUE = to be interpolated; logical "OR" is used to combine columns)
%  Tr           : Repetition time in seconds (default = 1 [sec])
%  InterpMethod : interpolation methods (either 'Lomb' or methods of interp1.m, e.g. 'linear', 'nearest', 'spline', 'cubic', etc.)
%  Prefix       : default = 'l'


%% Parse inputs

Img = cellstr(char(Img));
[pth,nme,ext] = shiFileParts(Img);

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'l';
end
if ~exist('Flag','var') || isempty(Flag)
    Flag = [];
end

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

if ~exist('Tr','var') || isempty(Tr)
    Tr = 1;
end
if ~exist('InterpMethod','var') || isempty(InterpMethod)
    InterpMethod = 'Lomb';
end

if ~exist('Mask','var') || isempty(Mask)
    applyMask = false;
else
    applyMask = true;
end

% default Lomb interpolation parameters (if InterpMethod == 'Lomb')

Oversampling_Frequency = [];
Highest_Frequency = [];
szChk = 500; % choose chunck size wisely to avoid memory error


%% get time stamps

TimeAll = (1:length(Img)) * Tr;

indSpike = false(size(Img));
if ischar(Flag)
    Flag = cellstr(char(Flag));
    for i = 1:length(Flag)
        indSpike = any([indSpike,readmatrix(Flag{i})>0]);
    end
elseif isnumeric(Flag)
    indSpike = any(Flag>0,2);
elseif islogical(Flag)
    indSpike = any(Flag,2);
elseif iscell(Flag)
    for i = 1:length(Flag)
        if ischar(Flag{i})
            indSpike = any([indSpike,readmatrix(Flag{i})>0],2);
        elseif isnumeric(Flag{i})
            indSpike = any([indSpike,Flag{i}>0],2);
        elseif islogical(Flag{i})
            indSpike = any([indSpike,Flag{i}],2);
        end
    end
else
    error('unknown format of Flag');
end

TimeX = TimeAll(~indSpike);
TimeY = TimeAll(indSpike);


%% main

V = spm_vol(char(Img));
All = spm_read_vols(V);
sz = size(All);

if applyMask
    [~,msk] = shiSpmMaskRead(Img{1},Mask,'incl');
else
    msk = true(sz(1:3));
end

All_reshape = reshape(All,prod(sz(1:3)),sz(4))'; % nVolume * nVoxel
msk_reshape = reshape(msk,numel(msk),1); % nVoxel * 1

All_msk = All_reshape(:,msk_reshape); % nVolume * nVoxelMasked
X_msk = All_msk(~indSpike,:); % nVolumeIncluded * nVoxelMasked

fprintf('Interpolate:');

if matches(lower(InterpMethod),'lomb')
    Y_msk = sinterp_lomb(TimeX, X_msk, TimeY, Oversampling_Frequency, Highest_Frequency, szChk);
else
    Y_msk = interp1(TimeX, X_msk, TimeY, InterpMethod);
end

All_msk(indSpike, :) = Y_msk; % nVolumeExcluded * nVoxelMasked

All_reshape(:,msk_reshape) = All_msk; % nVolume * nVoxelMasked
All_reshape(:,~msk_reshape) = 0; % nVolume * nVoxelOutsideMask

All = reshape(All_reshape',sz);


%% writing

fprintf(' writing...');

for i = 1:length(Img)
    Vo = V(i);
    Vo.fname = outImg{i};
    Vo.descrip = [Vo.descrip, sprintf('; Interpolate, method=%s', InterpMethod)];
    Vo.dt(1) = 16;
    spm_write_vol(Vo,single(All(:,:,:,i)));
end

fprintf(' Done.\n');



function Y = sinterp_lomb(TimeX,X,TimeY, oFreq, hiFreq, szChk)

nVox = size(X,2);

Y = nan(numel(TimeY),nVox);

p = 1;
q = szChk;

fprintf(' (Lomb)   0.00%%');

while p <= nVox
    q = min(q,nVox);
    Y(:,p:q) = shiStatLombInterp(TimeX, X(:,p:q), TimeY, oFreq, hiFreq);
    p = p + szChk;
    q = q + szChk;
    fprintf('\b\b\b\b\b\b\b%6.2f%%',min(100,q/nVox*100));
end