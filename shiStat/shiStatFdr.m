function [BH_corrected_p,BHY_corrected_p] = shiStatFdr(p,ignoreLowerTriag,ignoreDiagnal)
 
% returns the FDR correction results of raw p values using Benjamini-Hochberg-Yekutieli procedure (see wikipedia)
%
%       Benjamini, Yoav; Hochberg, Yosef (1995). "Controlling the false
%          discovery rate: a practical and powerful approach to multiple
%          testing". Journal of the Royal Statistical Society, Series B. 57
%          (1): 289-300.
%       Benjamini, Yoav; Yekutieli, Daniel (2001). "The control of the
%          false discovery rate in multiple testing under dependency".
%          Annals of Statistics. 29 (4): 1165-1188.
% 
%   
%   pval                  - scalar, vector, or matrix of raw p values
%   BH_corrected_p        - corrected p (B&H 1995)
%   BHY_corrected_p       - corrected p (B&Y 2001) (If the tests are under arbitrary dependence)
%
%    ###########
% by Zhenhao Shi @ 2016-12-21
%    ###########

if ~exist('ignoreLowerTriag','var') || isempty(ignoreLowerTriag)
    ignoreLowerTriag = false;
end
if ~exist('ignoreDiagnal','var') || isempty(ignoreDiagnal)
    ignoreDiagnal = false;
end
if size(p,1)==size(p,2) && ismatrix(p)
    idxNaN = false(size(p));
    if ignoreLowerTriag
        idxNaN = idxNaN | tril(ones(size(p)),-1)>0;
    end
    if ignoreDiagnal
        idxNaN = idxNaN | eye(size(p,1))>0;
    end
    p(idxNaN) = NaN;
end

[p_sorted,ix] = sort(p(:),'descend');
[~,ix_ix] = sort(ix);

n_Nan = sum(isnan(p(:)));

BH_corrected_p = nan(size(p));
BHY_corrected_p = nan(size(p));

FDR_BH_sorted = BH_corrected_p(:);
FDR_BHY_sorted = BHY_corrected_p(:);

N = numel(p)-n_Nan;
C = sum(1./(1:N));

FDR_BH_sorted(n_Nan+1) = p_sorted(n_Nan+1);
FDR_BHY_sorted(n_Nan+1) = min(1,p_sorted(n_Nan+1)).*C;

for k = n_Nan+2:n_Nan+N
    FDR_BH_sorted(k) = min(FDR_BH_sorted(k-1),p_sorted(k)*N/(N+1-k+n_Nan));
    FDR_BHY_sorted(k) = min(FDR_BHY_sorted(k-1),p_sorted(k)*N*C/(N+1-k+n_Nan));
end

BH_corrected_p(:) = FDR_BH_sorted(ix_ix);
BHY_corrected_p(:) = FDR_BHY_sorted(ix_ix);

