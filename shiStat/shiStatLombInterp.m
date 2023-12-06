function Y = shiStatLombInterp(TimeX, X, TimeY, OverFreq, HighFreq)

% performs interpolation based on the Lomb-Scargle periodogram
%
% TimeX     : available time points excluding flagged ones
% X         : data at available time points (one timepoint per row, one variable per column)
% TimeAll   : output timepoints (including flagged ones)
% OverFreq  : oversampling frequency (generally >= 4, default = 8)
% HighFreq  : allowing HighFreq*Nyquist_Limit as highest frequency (default = 1)
% 
% Adapted from Anish Mitra, October 25 2012

if ~exist('OverFreq','var') || isempty(OverFreq)
    OverFreq = 8;
end
if ~exist('HighFreq','var') || isempty(HighFreq)
    HighFreq = 1;
end

TimeX = TimeX(:);
TimeY = TimeY(:);

MEAN = mean(X,1);
X = bsxfun(@minus,X,MEAN);

[Nt, Nv] = size(X); % Number of time points x Number of voxels
T = max(TimeX) - min(TimeX); % Total time span

% calculate sampling frequencies
f = (1/(T*OverFreq) : 1/(T*OverFreq) : HighFreq*Nt/(2*T)).';

% angular frequencies and constant offsets
w = 2*pi*f;
tau = atan2(sum(sin(2*w*TimeX.'), 2), sum(cos(2*w*TimeX.'), 2))./(2*w);

% spectral power sin and cosine terms
cterm = cos(w*TimeX.' - repmat(w.*tau, 1, length(TimeX)));
sterm = sin(w*TimeX.' - repmat(w.*tau, 1, length(TimeX)));

D = reshape(X, 1, Nt, Nv);

cMult = bsxfun(@times, cterm, D);
C_final = bsxfun(@rdivide, sum(cMult, 2), sum(bsxfun(@power, cterm, 2), 2));

sMult = bsxfun(@times, sterm, D);
S_final = bsxfun(@rdivide, sum(sMult, 2), sum(bsxfun(@power, sterm, 2), 2));

% The inverse function to re-construct the original time series
T_rep = repmat(TimeY', [size(f, 1), 1, Nv]);
prod = bsxfun(@times, T_rep, w);
sw_p = bsxfun(@times, sin(prod), S_final);
cw_p = bsxfun(@times, cos(prod), C_final);
Y = reshape(sum(cw_p) + sum(sw_p), size(TimeY, 1), Nv);

% Normalize the reconstructed spectrum, needed when ofac > 1
norm_fac = std(Y)./std(X);
Y = bsxfun(@rdivide, Y, norm_fac);
Y = bsxfun(@plus,Y,MEAN);