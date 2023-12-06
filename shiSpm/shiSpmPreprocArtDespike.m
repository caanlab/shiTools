function outImg = shiSpmPreprocArtDespike(Img,WinSize,Thres,Prefix,existAction)

% INPUTS
%    Img:      a list of images, in a single session.
%    WinSize:  must be odd number; default = 17
%    Thres:    clip threshold used to despike the data, in units of percentage signal change relative to the mean image; default = 4 (i.e. 4% deviation from mean)
%
%    FiltType has 4 options, only using option 3 here:
%       3. No high pass filtering is done. Despiking is done based
%          on a WinSize-point moving average of unfiltered data.
%
% adapted from ArtRepair art_despike
%


Img = cellstr(char(Img));
[pth,nme,ext] = shiFileParts(Img);

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'k';
end
if ~exist('WinSize','var') || isempty(WinSize)
    WinSize = 17;
end
if ~exist('Thres','var') || isempty(Thres)
    Thres = 4;
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

if mod(WinSize,2) == 0   % check that filter length is odd.
    disp('Warning from art_despike: Filter length must be odd')
    return
end


%%

nImg = size(Img,1);

V = spm_vol(char(Img));
Y = spm_read_vols(V);
Ymean = mean(Y,4);



%%

Lag = (WinSize-1)/2;  % filtered value at 9 corresponds to a peak at 5.
% Convert despike threshold in percent to fractional mulitplier limits
SpkLimitUpper = 1 + 0.01*Thres;
SpkLimitLower = 1 - 0.01*Thres;


% FILTER AND DESPIKE IMAGES
% Process all the scans sequentially
% Start and End are padded by reflection, e.g. sets data(-1) = data(3).
% Initialize lagged values for filtering with reflected values
% Near the end, create forward values for filtering with reflected values.

% Find mean image


% Initialize arrays with reflected values.
Y4 = zeros(WinSize,size(Ymean,1),size(Ymean,2),size(Ymean,3));
for i = 1:(Lag+1)
    i2 = i + Lag;
    Y4(i2,:,:,:) = (Y(:,:,:,i));
    i3 = Lag + 2 - i;
    if i > 1   % i=1 then i3 = i2.
        Y4(i3,:,:,:) = Y4(i2,:,:,:);
    end
end
%  Start up clipping is done here
movmean = squeeze(mean(Y4,1));
for i = 1:WinSize
   Y4s = squeeze(Y4(i,:,:,:));
   Y4s = min(Y4s,SpkLimitUpper*movmean);
   Y4s = max(Y4s,SpkLimitLower*movmean);
   Y4(i,:,:,:) = Y4s;
end

% Main Loop
% Speed Note: Use Y4(1,:,:,:) = spm_read_vols(P(1));  % rows vary fastest
for i = (WinSize+1)/2 : nImg+(WinSize-1)/2
    if i <= nImg
        Y4(WinSize,:,:,:) = (Y(:,:,:,i));
    else   % Must pad the end data with reflected values.
        Y4(WinSize,:,:,:) = (Y(:,:,:,2*nImg-i));
    end
    %  Incremental clipping is done here
    movmean = squeeze(mean(Y4,1));
    % This lag is from FiltType = 3
    Y4s = squeeze(Y4(WinSize-Lag,:,:,:));  % centered for despike only
    Y4s = min(Y4s,SpkLimitUpper*movmean);
    Y4s = max(Y4s,SpkLimitLower*movmean);
    Yn2 = squeeze(Y4s);

    % Prepare the header for the filtered volume, with lag removed.
    Vo = V(1);
    Vo.fname = outImg{i-Lag};
    Vo.descrip = [Vo.descrip, sprintf('; art_despike, WinSize=%d, Thres=%d', WinSize, Thres)];
    spm_write_vol(Vo,Yn2); 

    for js = 1:WinSize-1
        Y4(js,:,:,:) = Y4(js+1,:,:,:);
    end 
end
