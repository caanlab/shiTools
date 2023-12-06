function [z,mu,sigma] = nanzscore(x,flag,dim)
% returns standardized z scores ignoring NaNs
%   Z = nanzscore(X) returns a centered, scaled version of X, the same size as X.
%   For vector input X, Z is the vector of z-scores (X-nanmean(X)) ./ nanstd(X). For
%   matrix X, z-scores are computed using the nanmean and standard deviation
%   along each column of X.  For higher-dimensional arrays, z-scores are
%   computed using the nanmean and standard deviation along the first
%   non-singleton dimension.
%
%   The columns of Z have sample nanmean zero and sample standard deviation one
%   (unless a column of X is constant, in which case that column of Z is
%   constant at 0).
%
%   [Z,MU,SIGMA] = nanzscore(X) also returns nanmean(X) in MU and nanstd(X) in SIGMA.
%
%   [...] = nanzscore(X,1) normalizes X using nanstd(X,1), i.e., by computing the
%   standard deviation(s) using N rather than N-1, where N is the length of
%   the dimension along which nanzscore works.  nanzscore(X,0) is the same as
%   nanzscore(X).
%
%   [...] = nanzscore(X,FLAG,DIM) standardizes X by working along the dimension
%   DIM of X. Pass in FLAG==0 to use the default normalization by N-1, or 1
%   to use N.
%
%   See also nanmean, nanstd.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2010/03/16 00:18:33 $

% [] is a special case for nanstd and nanmean, just handle it out here.
if isequal(x,[]), z = []; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's nanmean and sd, and standardize it
mu = nanmean(x,dim);
sigma = nanstd(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);

