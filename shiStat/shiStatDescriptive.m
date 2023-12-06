function stat = shiStatDescriptive(X)

% returns descriptive statistics summary for each column of X
% 
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########

stat = struct([]);
for i = 1:size(X,2)
    x=X(:,i);
    stat(i).N = length(x);
    stat(i).N_valid = sum(~isnan(x));
    stat(i).N_missing = sum(isnan(x));
    stat(i).Mean = nanmean(x);
    stat(i).Median = nanmedian(x);
    stat(i).SD = nanstd(x);
    stat(i).SEM = nanstd(x)/sqrt(sum(~isnan(x)));
    stat(i).Min = min(x);
    stat(i).Max = max(x);
    [~,temp_p,~,temp_stat] = ttest(x);
    stat(i).T = temp_stat.tstat;
    stat(i).Df = temp_stat.df;
    stat(i).P = temp_p;
    [~,temp_p,temp_stat] = swtest(x);
    stat(i).Normality_SWtest_W = temp_stat;
    stat(i).Normality_SWtest_P = temp_p;
end