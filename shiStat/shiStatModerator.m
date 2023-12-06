function [Int_Beta,Int_t,Int_df,Int_p,stats] = shiStatModerator(y,X,M,C,verbose)
 
% conducts moderator analysis
%
% [Int_Beta,Int_t,Int_df,Int_p,stats] = shiStatModerator(y,X,M)
%
%   y           - n-by-1 vector of dependent variable
%   X           - n-by-1 vector of independent variable
%   M           - n-by-1 vector of moderator
%   Int_Beta    - beta value of interaction term
%   Int_t       - t statistic value of interaction term
%   Int_df      - degree of freedom for t statistic
%   Int_p       - significance of interaction term
%   stats       - other statistics (R^2, Radj^2, ts, fs)
% 
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########

if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

if ~exist('C','var') || isempty(C)
    C = [];
else
    if size(C,1) ~= size(y,1)
        error('number of observations not matched');
    end
end

if size(y,2) ~= 1
    error('y is not a vector');
end

if ~isequal(size(y),size(M)) || ~isequal(size(M),size(X))
    error('number of observations not matched');
end

y = nanzscore(y);
M = nanzscore(M);
X = nanzscore(X);
Int = nanzscore(M.*X);

if size([M,X,Int,C],1) <= size([M,X,Int,C],2)
    warning('no enough observations');
    [Int_Beta,Int_t,Int_df,Int_p,stats] = deal(NaN,NaN,NaN,NaN,struct([]));
    return;
end

stats = regstats(y,[M,X,Int,C],'linear',{'rsquare','adjrsquare','tstat','fstat'});

Int_Beta = stats.tstat.beta(4);
Int_t = stats.tstat.t(4);
Int_df = stats.tstat.dfe;
Int_p = stats.tstat.pval(4);

if verbose
    fprintf('\n  beta_int=%+.2f, N=%d, t(%d)=%+.2f, p=%g\n\n', Int_Beta, sum(~isnan(sum([y,M,X,Int,C],2))),Int_df, Int_t, Int_p);
end
