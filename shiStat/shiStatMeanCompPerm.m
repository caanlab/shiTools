function p = shiStatMeanCompPerm(Vec1,Vec2,n_iter)

% compares the mean of two independent samples using permutation test
%
% zhenhao shi

Vec1 = Vec1(~isnan(Vec1));
Vec2 = Vec2(~isnan(Vec2));

n1 = numel(Vec1);
n2 = numel(Vec2);

Vec_All = [Vec1(:);Vec2(:)];

m_diff = mean(Vec1(:)) - mean(Vec2(:));

m_diff_perm = nan(n_iter,1);

for i = 1:n_iter
    
    Vec_perm = Vec_All(randperm(n1+n2));
    
    m_diff_perm(i) = mean(Vec_perm(1:n1)) - mean(Vec_perm(n1+1:end));
    
end

p = sum( abs(m_diff) < abs(m_diff_perm) ) / n_iter;
