function QA = shiSpmPreprocQa(IMG,MASK)

% performs quality assurance
%
% see shiSpmPreprocQa.pptx
%

if ~exist('MASK','var') || isempty(MASK)
    MASK = 1/8;
end

IMG = cellstr(char(IMG));
[pth,nme,ext] = fileparts(IMG{1});
p = length(IMG);
V = spm_vol(char(IMG));


cwd=pwd;
if ~isempty(fileparts(V(1).fname))
    cd(fileparts(V(1).fname));
end

if isnumeric(MASK)
    Y = spm_read_vols(V);
    Ym = isfinite(mean(Y,4)) & mean(Y,4)>mean(Y(:))*MASK;
    Vm = V(1);
    Vm.fname = fullfile(pth,['QaMask_',nme,ext]);
    Vm.dt = [2 0];
else
    MASK = cellstr(char(MASK));
    [Vm,Ym,~,Y] = shiSpmMaskRead(IMG,MASK);
    Ym = isfinite(mean(Y,4)) & Ym;
    Vm.fname = fullfile(pth,['QaMask_',nme,ext]);
end

spm_write_vol(Vm,Ym);
MASK = Vm.fname;
Y = reshape(Y,prod(V(1).dim(1:3)),p)';
indY = find(Ym);
Y = Y(:,indY); 


%

disp('QA: Step 1 ...');
Y01 = fDetrend(Y);
Y02 = fTmMean(Y);
Y03 = fOdEvDiff(Y);
Y04 = fSpMean(Y);

disp('QA: Step 2 ...');
Y05 = fTmZscore(Y01);
Y06 = fTrendFit(Y,Y01);
Y07 = fTmSd(Y01);
Y08 = fSpMean(Y02);
Y09 = fSpSd(Y03);
Y10 = fDetrend(Y04);

disp('QA: Step 3 ...');
[imgRPV,paramFWHM] = fSmoothEst(V,Y05,MASK,indY,'QaRpv_');
Y12 = fDrift(Y06);
Y13 = fTmSnr(Y02,Y07);
Y14 = fFluct(Y02,Y07);
Y15 = fSpSnr(Y08,Y09,p);
Y16 = fTmSd(Y10);
Y17 = fTrendFit(Y04,Y10);

disp('QA: Step 4 ...');
imgDRIFT = fImgWrite(V,Y12,indY,'QaDrift_');
imgTEMPSNR = fImgWrite(V,Y13,indY,'QaTempSnr_');
Y18 = fSpMean(Y13);
imgFLUCT = fImgWrite(V,Y14,indY,'QaFluct_');
Y19 = fFluct(Y08,Y16);
paramSPATSNR = Y15;
Y20 = fDrift(Y17);

disp('QA: Step 5 ...');
Y21 = fRms(paramFWHM);
paramTEMPSNR = Y18;
paramFLUCT = Y19;
paramDRIFT = Y20;

disp('QA: Step 6 ...');
paramFWHMSUMM = Y21;

disp('QA: done.');
QA = struct(...
    'RPV_img',      imgRPV, ... % resels per voxel (image)
    'FWHM_x',       paramFWHM(1), ... % smoothness on x, in mm
    'FWHM_y',       paramFWHM(2), ... % smoothness on y, in mm
    'FWHM_z',       paramFWHM(3), ... % smoothness on z, in mm
    'FWHM_rms',     paramFWHMSUMM, ...  % overall smoothness, in mm
    'Drift_img',    imgDRIFT, ... % percentage signal drift (image)
    'Drift',        paramDRIFT, ... % percent signal drift
    'TempSnr_img',  imgTEMPSNR, ... % tSNR (image)
    'TempSnr',      paramTEMPSNR, ... % tSNR
    'SpatSnr',      paramSPATSNR, ... % spatial SNR
    'Fluct_img',    imgFLUCT, ... % percent signal fluctuation (image)
    'Fluct',        paramFLUCT ... % percent signal fluctuation
    );

save(['QA_',nme,'.mat'],'QA');
cd(cwd);

%
function Y = fDetrend(Y)
Y = spm_detrend(Y,2);

function Y = fTrendFit(Y_Raw,Y_Detrend)
Y = Y_Raw - Y_Detrend;

function Y = fTmMean(Y)
Y = mean(Y,1);

function Y = fTmSd(Y)
Y = std(Y,[],1);

function Y = fTmZscore(Y)
Y = zscore(Y,[],1);

function Y = fSpMean(Y)
Y = nanmean(Y,2);

function Y = fSpSd(Y)
Y = nanstd(Y,[],2);

function Y = fOdEvDiff(Y)
Y = sum(Y(1:2:2*floor(size(Y,1)/2)-1,:)-Y(2:2:2*floor(size(Y,1)/2),:),1);

function Y = fRms(Y)
Y = sqrt(mean(Y(:).^2));

function Y = fDrift(Y)
Y = range(Y,1)./mean(Y,1).*100;

function Y = fFluct(Mean,Sd)
Y = Sd./Mean.*100;

function Y = fTmSnr(Mean,Sd)
Y = Mean./Sd;

function Y = fSpSnr(Mean,Sd,p)
Y = Mean./(Sd./sqrt(p));

function [OutName,paramFWHM] = fSmoothEst(V,Y,MASK,indY,Prefix)
RAND = [shiTime,sprintf('_%d_',round(rand()*10000))];
VoxSize = [sqrt(sum(V(1).mat(1:3,1).^2)),sqrt(sum(V(1).mat(1:3,2).^2)),sqrt(sum(V(1).mat(1:3,3).^2))];
[pth,nme,ext] = fileparts(V(1).fname);
OutName = fullfile(pth,[Prefix,nme,ext]);
Yo = nan(V(1).dim);
ZNAME = cell(size(Y,1),1);
for i = 1:size(Y,1)
    Yo(indY) = Y(i,:);
    Vo = V(1);
    Vo.fname = fullfile(pth,['_tempQa_',RAND,nme,num2str(i,'%03d'),ext]);
    ZNAME{i} = Vo.fname;
    spm_write_vol(Vo,Yo);
end
[paramFWHM,VRpv] = spm_est_smoothness(char(ZNAME),MASK);
paramFWHM = paramFWHM(:).*VoxSize(:);
movefile(VRpv.fname,OutName);
for i = 1:size(Y,1)
    delete(ZNAME{i});
end
if ~strcmpi(ext,'.img')
    [zpth,znme] = shiFileParts(ZNAME);
    ZNAME2 = shiStrConcat(zpth,filesep,znme,'.hdr');
    for i = 1:size(Y,1)
        delete(ZNAME2{i});
    end
end

function OutName = fImgWrite(V,Y,indY,Prefix)
[pth,nme,ext] = fileparts(V(1).fname);
OutName = fullfile(pth,[Prefix,nme,ext]);
V = V(1);
Yo = nan(V.dim);
Yo(indY) = Y;
V.fname = OutName;
V.descrip = [V.descrip,'; QA'];
spm_write_vol(V,Yo);
