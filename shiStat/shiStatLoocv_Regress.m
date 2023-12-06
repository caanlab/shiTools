function [MSE,R,Y_fitted] = shiStatLoocv_Regress(Y,X,n_Iter)

% performs leave-one-out cross-validation and returns mean squared error of prediction, correlation between actual Y and predicted Y, and values of predicted Y
%
% zhenhao shi

if size(Y,2)~=1 || size(Y,1) ~= size(X,1)
    error('unmatched dimension of Y and X');
end

ListWise = ~any(isnan([Y,X]),2);
Y = Y(ListWise,1);
X = X(ListWise,:);

[MSE.MSE,R.R,Y_fitted] = shiStatLoocv_Regress_func(Y,X);

fprintf('\nLOOCV 00.00%%');
for j = 1:n_Iter
    [MSE.MSE_perm(j,1),R.R_perm(j,1)] = shiStatLoocv_Regress_func(Y(randperm(size(Y,1))),X);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\bLOOCV %05.2f%%',j/n_Iter*100);
end

MSE.p = sum(MSE.MSE_perm<MSE.MSE)/n_Iter;
R.p = sum(R.R_perm>R.R)/n_Iter;

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bLOOCV DONE\n');


function [MSE,R,Y_fitted] = shiStatLoocv_Regress_func(Y,X)

Y_fitted = nan(size(Y));

for k = 1:size(X,1)
    Ioo = setdiff(1:size(X,1),k);
    Yoo = Y(Ioo,1);
    Xoo = X(Ioo,:);
    beta = regstats(Yoo,Xoo,'linear',{'beta'});
    beta = beta.beta;
    Y_fitted(k) = [1,X(k,:)]*beta;
end

MSE = mean((Y-Y_fitted(:)).^2);
R = corr(Y,Y_fitted(:));

