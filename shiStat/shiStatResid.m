function resid_Y = shiStatResid(Y,Cov,addConstantToX,addMeanToResid)

% regresses covariate(s) out of a variable and returns residual
% 
% resid_Y = shiStatResid(Y,Cov)
% resid_Y = shiStatResid(Y,Cov,addConstantFirst)
%   returns the residuals of Y after regressing out Cov
% 
%   Y    - vector, input dependent variable
%   Cov  - matrix of covariates, each column is a covariate
%   addConstantToX  - whether to add constant term to the first column
%                       (default = true)
% 
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########

if isempty(Cov)
    resid_Y = Y;
    return;
end

if nargin < 3 || isempty(addConstantToX)
    addConstantToX = true;
end
if nargin < 4 || isempty(addMeanToResid)
    addMeanToResid = true;
end

if size(Y,1) ~= size(Cov,1)
    error('number of observation not matched');
end

resid_Y = nan(size(Y));

if addConstantToX || ~isequal(ones(size(Y,1),1),Cov(:,1))
    Cov = [ones(size(Y,1),1),Cov];
end

for i = 1:size(Y,2)
    [~,~,resid_Y(:,i)] = regress(Y(:,i),Cov);
    if addMeanToResid
        resid_Y(:,i) = resid_Y(:,i)+nanmean(Y(:,i));
    end
end

