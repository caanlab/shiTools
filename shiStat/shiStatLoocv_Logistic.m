function [ROC,Prob_fitted] = shiStatLoocv_Logistic(Y,X,n_Iter)

% performs leave-one-out cross-validation for logistic regression and returns ROC results and fitted probabilities
%
% zhenhao shi

if size(Y,2)~=1 || size(Y,1) ~= size(X,1)
    error('unmatched dimension of Y and X');
end

if numel(unique(Y(:))) ~= 2
    error('Y must be 0/1');
elseif islogical(Y)
    Y = Y.*1;
elseif min(Y(:)) ~= 0 || max(Y(:)) ~= 1
    error('Y must be 0/1');
end



ListWise = ~any(isnan([Y,X]),2);
Y = Y(ListWise,1);
X = X(ListWise,:);

[ROC.Coord_x,ROC.Coord_y,ROC.AUC,Prob_fitted,ROC.MaxAcc] = shiStatLoocv_Logistic_func(Y,X);

fprintf('\nLOOCV 00.00%%');
for j = 1:n_Iter
    [ROC_Coord_x_perm{j,1},ROC_Coord_y_perm{j,1},ROC_AUC_perm(j,1),~,ROC_MaxAcc_perm(j,1)] = shiStatLoocv_Logistic_func(Y(randperm(length(Y))),X);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\bLOOCV %05.2f%%',j/n_Iter*100);
end
ROC.Coord_x_perm = ROC_Coord_x_perm;
ROC.Coord_y_perm = ROC_Coord_y_perm;
ROC.AUC_perm = ROC_AUC_perm;
ROC.MaxAcc_perm = ROC_MaxAcc_perm;

if n_Iter > 0
    ROC.pval_AUC = sum(ROC.AUC<ROC.AUC_perm)/n_Iter;
    ROC.pval_MaxAcc = sum(ROC.MaxAcc<ROC.MaxAcc_perm)/n_Iter;
end

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bLOOCV DONE\n');


function [ROC_x,ROC_y,ROC_AUC,Score_fitted,MaxAcc] = shiStatLoocv_Logistic_func(Y,X)

Score_fitted = nan(size(Y));

for k = 1:size(X,1)
    
    Ioo = setdiff(1:size(X,1),k);
    Yoo = Y(Ioo,1);
    Xoo = X(Ioo,:);

    MDLoo = fitglm(Xoo,Yoo,'linear','Distribution','binomial');
    
    Xnew = num2cell(X(k,:));
    Score_fitted(k) = feval(MDLoo,Xnew{:});

end

[ROC_x,ROC_y,thres,ROC_AUC,optpnt] = perfcurve(Y,Score_fitted,1);
thres = thres(ROC_x==optpnt(1) & ROC_y==optpnt(2));
MaxAcc = sum((Y>.5) == (Score_fitted>=thres(end)))/length(Y);
ROC_x=ROC_x';
ROC_y=ROC_y';


