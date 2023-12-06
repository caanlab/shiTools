function [PseudoF,Df1,Df2,Pval,PseudoF_Perm] = shiStatMdmr(DistMat,X_Iv,X_Cov,n_Perm)

% performs multivariate distance-based matrix regression and returns pseudo-F value, degrees of freedom, and p value from permutation test
%
% zhenhao shi

if nargin < 4 || isempty(n_Perm)
    n_Perm = 2000;
end
if nargin < 3 || isempty(X_Cov)
    X_Cov = [];
end
if n_Perm <= 0
    error('number of permutation must be greater than 0');
end

if size(DistMat,1) ~= size([X_Iv,X_Cov],1)
    error('sample size not matched');
elseif isequal(DistMat,DistMat')
    error('distance matrix not symetric');
end

N = size(DistMat,1);

if isempty(X_Cov) || ~isequal(X_Cov(:,1),ones(N,1))
    X_Cov = [ones(N,1),X_Cov];
end

X0 = X_Cov(:,1);
X1 = zscore(X_Cov(:,2:end));
X2 = zscore(X_Iv);

if rank([X0,X1,X2])<size([X0,X1,X2],2)
    error('X is rank deficient');
end

PseudoF = shi_func_PseudoF(DistMat,X0,X1,X2);
Df1 = size(X2,2);
Df2 = N-size([X0,X1,X2],2);

beta_X2 = nan(size([X0,X1],2),size(X2,2));
resid_X2 = nan(size(X2));
fitted_X2 = nan(size(X2));
for i = 1:size(X2,2)
    [beta_X2(:,i),~,resid_X2(:,i)] = regress(X2(:,i),[X0,X1]);
    fitted_X2(:,i) = [X0,X1]*beta_X2(:,i);
end

PseudoF_Perm = nan(n_Perm,1);
for p = 1:n_Perm
    perm_X2 = fitted_X2 + resid_X2(randperm(N),:);
    PseudoF_Perm(p,1) = shi_func_PseudoF(DistMat,X0,X1,perm_X2);
end

Pval = sum(PseudoF<PseudoF_Perm)/n_Perm;


function PseudoF = shi_func_PseudoF(DistMat,X0,X1,X2)

G = shi_func_D2G(DistMat);
H = shi_func_X2H([X0,X1,X2]);
H2 = H - shi_func_X2H([X0,X1]);
xx1 = trace(H2*G) / (size(X2,2));
xx2 = trace((eye(size(X0,1))-H)*G) /(size(X0,1)-size([X0,X1,X2],2));
PseudoF = xx1/xx2;

function G = shi_func_D2G(D)
N = size(D,1);
ONE = ones(N,1);
EYE = eye(N);
A = -D.^2/2;
G = (EYE - ONE*ONE'./N) * A * (EYE - ONE*ONE'./N);

function H = shi_func_X2H(X)
H = X * (X'*X)^(-1) * X';