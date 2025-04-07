function outImg = shiSpmPreproc3dDespike(Img,Mask,c1,c2,cOrder,Prefix,existAction)
%
% performs AFNI's 3dDespike (using the NEW method)
%
% Img      : Raw images
% MaskImg  : Inclusive mask image
% c1,c2    : Cut range for despiking (default: 2.5 and 4.0)
% cOrder   : Curve fit order. Sets the number of sin and cosine terms to fit. (Default is number of volumes / 30).
% Prefix   : defalt = 'k'
%
% 2024.02.26


sz_chk = 300;


%% Parse inputs

Img = cellstr(char(Img));
[pth,nme,ext] = shiFileParts(Img);
if isempty(pth{1}), pth = pwd; end

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'k';
end
if ~exist('c1','var') || isempty(c1)
    c1 = 2.5;
end
if ~exist('c2','var') || isempty(c2)
    c2 = 4.0;
end
if ~exist('cOrder','var') || isempty(cOrder)
    cOrder = round(length(Img)/30);
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

if ~exist('Mask','var') || isempty(Mask)
    applyMask = false;
else
    applyMask = true;
end


%%

V = spm_vol(char(Img));
Yin = spm_read_vols(V);

if applyMask
    [~,Y_Mask] = shiSpmMaskRead(Img{1},Mask,'incl');
else
    Y_Mask = true(size(Yin(:,:,:,1)));
end

PolyOrder = 2; %% as in AFNI 3ddespike

X = DsnMtrx(size(Yin,4), cOrder, PolyOrder);

sdY = std(Yin,[],4);


msk = Y_Mask(:)>0 & sdY(:)>0;
Y = reshape(Yin,prod(size(Yin,1:3)),size(Yin,4))';
Y = Y(:,msk);

n_chk = ceil(size(Y,2)/sz_chk);

i1 = 1;
i2 = min(size(Y,2),sz_chk);


cnt = 0;
despiked = 0;
fprintf('\n3dDespike:       ');

invX = (X.' * X) \ X.';

while i1 <= size(Y,2)

    cnt = cnt+1;

    xY = Y(:,i1:i2);
    % xB = shiStatLeastAbsDevReg(xY, X); % LAD regression using IRLS
    xB = DES_solve(xY,invX); % 3ddespike NEW algorithm

    % residual from predicted
    xY_predict = X * xB;
    xY_resid = xY - xY_predict;      % Get actual value minus predicted values

    % scaled median absolute deviation
    xY_sigma = sqrt(pi/2) .* median(abs(xY_resid));

    % standardized residual
    xY_sresid = xY_resid ./ xY_sigma;

    % identify and fix spikes
    ind_spike = abs(xY_sresid) > c1;
    xY_sresid(ind_spike) = sign(xY_sresid(ind_spike)) .* (  c1 + (c2-c1) .* tanh( (abs(xY_sresid(ind_spike))-c1) ./ (c2-c1) )  );

    % reconstruct xY
    xY_despike = xY_sresid .* xY_sigma + xY_predict;
    xY_despike(:,any(isnan(xY_despike))) = xY(:,any(isnan(xY_despike)));

    Y(:,i1:i2) = xY_despike;
    despiked = despiked + sum(ind_spike(:));

    i1 = i1 + sz_chk;
    i2 = min(size(Y,2),i2 + sz_chk);

%     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%5.1f%%',cnt/n_chk*100);
    fprintf('\b\b\b\b\b\b%5.1f%%',cnt/n_chk*100);
    
end

pct_dspkd = despiked/sum(msk(:))/size(Y,1)*100;
fprintf('DONE: despiked %f5.2%% values\n',pct_dspkd)

Yout = zeros(prod(size(Yin,1:3)),size(Yin,4))';
Yout(:,msk) = Y;
Yout = reshape(Yout',size(Yin));

for i = 1:length(Img)
    Vo = V(i);
    Vo.fname = outImg{i};
    Vo.descrip = [Vo.descrip, sprintf('; 3ddespikeNEW, c1=%g, c2=%g, corder=%d, pct_dspkd=%g', c1, c2, cOrder,pct_dspkd)];
    spm_write_vol(Vo,Yout(:,:,:,i));
end


%%

function X = DsnMtrx(nImg, cOrder, polyOrder)

nVar = 2 * cOrder + polyOrder + 1;

X = nan(nImg, nVar);

Poly = 0:polyOrder;
x = (1:nImg)';
for i = 1:length(Poly)
    X(:,i) = x.^Poly(i);
end

cnt = length(Poly) + 1;
for i = 1:cOrder
    X(:, cnt) = sin(((1:nImg)' * 2 * pi * i) / nImg);
    cnt = cnt + 1;
    X(:, cnt) = cos(((1:nImg)' * 2 * pi * i) / nImg);
    cnt = cnt + 1;
end


%%

function [med,mad] = mead9(xY)

[ntim,nvox] = size(xY);
[med,mad] = deal(nan(ntim,nvox));

ja = [ones(1,4), 1:ntim-8, repmat(ntim-8,1,4)];
jb = ja + 8;

for jj = 1:ntim
    med(jj,:) = median(xY(ja(jj):jb(jj),:));
    mad(jj,:) = median(abs( xY(ja(jj):jb(jj),:) - med(jj,:) ));
end


%%

function xY = DES_despike9(xY)

[zme,zma] = mead9(xY);

ind = abs(xY - zme) > median(zma) * 6.789; % threshold value copied from 3ddespike

xY(ind) = zme(ind);

if any(ind(:))
    fprintf('♂');
end


%%

function coef = DES_solve(xY,invX)

coef = invX * DES_despike9(xY);

