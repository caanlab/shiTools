function [Beta, Pval, Stat, Boot] = shiStatRegressBoot(Y,X,n_iter)

% conducts bootstrap test of linear regression (may be revised or removed in future)
%
% zhenhao shi


nObs = size(X,1);
Const = ones(nObs,1);
if isequal(X(:,1),Const)
    X = X(:,2:end);
end
nVar = size(X,2);
X = nanzscore(X);
Y = nanzscore(Y);



stat = regstats(Y,X,'linear',{'tstat','fstat','adjrsquare'});
Beta = stat.tstat.beta;
tval = stat.tstat.t;
tstat = stat.tstat;
fstat = stat.fstat;
adjrsquare = stat.adjrsquare;
rho2 = nan(nVar+1,1);
for i = 1:nVar
    rho2(i+1,1) = partialcorr(Y,X(:,i),X(:,setdiff(1:nVar,i)));
end


beta_perm = nan(n_iter,nVar);
tval_perm = nan(n_iter,nVar);
rho2_perm = nan(n_iter,nVar);
pval_beta = nan(nVar+1,1);
pval_tval = nan(nVar+1,1);
pval_rho2 = nan(nVar+1,1);
fval_perm = nan(n_iter,1);
for i = 1:nVar
    X_reduced = X(:,setdiff(1:nVar,i)); % X
    X_Z = X(:,i);
    [beta_reduced,~,resid_reduced,~] = regress(Y,[Const,X_reduced]);   % [b_0 b] and R_(Y,X)
    Y_reduced_fitted = [Const,X_reduced]*beta_reduced;
    for j = 1:n_iter
        RandInd = ceil(rand(1,nObs)*nObs);
        resid_perm = resid_reduced(RandInd);
        Y_perm = Y_reduced_fitted + resid_perm;
        stat_perm = regstats(Y_perm,X,'linear',{'tstat'});
        beta_perm(j,i) = stat_perm.tstat.beta(i+1);
        tval_perm(j,i) = stat_perm.tstat.t(i+1);
        rho2_perm(j,i) = partialcorr(resid_perm,X_Z,X_reduced)^2;
    end
    pval_beta(i+1) = sum(Beta(i+1)<beta_perm(:,i))/n_iter;
    pval_tval(i+1) = sum(tval(i+1)<tval_perm(:,i))/n_iter;
    pval_rho2(i+1) = sum(rho2(i+1)<rho2_perm(:,i))/n_iter;
end

Pval = pval_beta;


for j = 1:n_iter
    RandInd = ceil(rand(1,nObs)*nObs);
    Y_perm = Y(RandInd);
    stat_perm = regstats(Y_perm,X,'linear',{'fstat'});
    fval_perm(j,1) = stat_perm.fstat.f;
end
pval_model_fval = sum(fstat.f < fval_perm)/n_iter;


Stat.tstat = tstat;
Stat.fstat = fstat;
Stat.partialcorr2 = rho2;
Stat.n_iter = n_iter;
Stat.pval_beta = pval_beta;
Stat.pval_tval = pval_tval;
Stat.pval_rho2 = pval_rho2;
Stat.pval_fval = pval_model_fval;
Stat.adjrsquare = adjrsquare;


Boot.beta_perm = beta_perm;
Boot.tval_perm = tval_perm;
Boot.rho2_perm = rho2_perm;
Boot.fval_perm = fval_perm;
    
        



