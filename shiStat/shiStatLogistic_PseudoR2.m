function [R2_Nagelkerke,R2_CoxSnell,R2_McFadden,R2_McFadden_adj] = shiStatLogistic_PseudoR2(Y,X)

% returns pseudo R2 for logistic regression
%
% Zhenhao Shi

ALL = [X,Y];
ALL = shi_deNaN(ALL);
[X,Y] = deal(ALL(:,1:end-1),ALL(:,end));

n_sub = length(Y);
n_var = size(X,2);

if exist('fitglm','file') || exist('fitglm','builtin')

    mdl = fitglm(X,Y,'linear','distribution','binomial');
    warning off;
    mdl0 = fitglm(zeros(size(Y)),Y,'linear','distribution','binomial');
    warning on;

    loglik = mdl.LogLikelihood;
    loglik0 = mdl0.LogLikelihood;

else

    b = glmfit(X,Y,'binomial','link','logit');
    lin = [ones(size(Y)), X] * b;
    loglik = sum(Y .* lin - log(1 + exp(lin)));

    warning off;
    b0 = glmfit(zeros(size(Y)),Y,'binomial','link','logit');
    loglik0 = sum(Y .* b0(1) - log(1 + exp(b0(1))));
    warning on;

end

% Nagelkerke R^2
R2_Nagelkerke = (1 - (exp(loglik0)/exp(loglik))^(2/n_sub) ) / (1 - exp(loglik0)^(2/n_sub));

% Cox and Snell R^2
R2_CoxSnell = 1 - (exp(loglik0)/exp(loglik))^(2/n_sub);

% McFadden R^2 
R2_McFadden = 1 - loglik/loglik0; % equivalent to 1 - dev/dev0
R2_McFadden_adj = 1 - (loglik - n_var)/loglik0; % counterpart to R2adj

