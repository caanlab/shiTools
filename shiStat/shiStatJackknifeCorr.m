function [DynCorr_raw, DynCorr_z] = shiStatJackknifeCorr(x,y,WinSize,CorrType)

% performes jackknife and d-jackknife dynamic correlation
% Thompson, W. H., Richter, C. G., Plavén-Sigray, P., & Fransson, P. (2018). Simulations to benchmark time-varying connectivity methods for fMRI. PLoS computational biology, 14(5), e1006196.
% Xie, H., Zheng, C. Y., Handwerker, D. A., Bandettini, P. A., Calhoun, V. D., Mitra, S., & Gonzalez-Castillo, J. (2019). Efficacy of different dynamic functional connectivity methods to capture cognitively relevant information. Neuroimage, 188, 502-514.

if nargin < 2
    y=x;
end

if ~exist('WinSize','var') || isempty(WinSize) || ~isfinite(WinSize) || WinSize<0
    WinSize = 1;
end
if ~exist('CorrType','var') || isempty(CorrType)
    CorrType = 'pearson';
end

if mod(WinSize,2) ~= 1
    error('WinSize must be an odd number');
end

DynCorr_raw = nan(size(x,2),size(y,2),size(x,1));

Win = (1:WinSize)-(1+WinSize)/2;

for i = 1:size(x,1)
    ind = setdiff(1:size(x,1),Win+i);
    DynCorr_raw(:,:,i) = corr(x(ind,:),y(ind,:),'rows','pairwise','type',CorrType);
end

DynCorr_z = zscore(DynCorr_raw,[],3);