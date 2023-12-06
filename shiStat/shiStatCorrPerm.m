function [r_real, pval, r_perm] = shiStatCorrPerm(Y,X,n_iter)

% performs permutation test of Pearson correlation
%

r_perm = nan(n_iter,1);
Y=Y(:);
X=X(:);

for i = 1:n_iter
    r_perm(i) = corr(Y,X(randperm(length(X))),'rows','pairwise');
end

r_real = corr(Y,X,'rows','pairwise');

pval = ( sum(abs(r_real)<r_perm) + sum(-abs(r_real)>r_perm) ) / n_iter;


