function [PV,thr]=shiStatPredictiveValue(z,x,n_Point,do_Reverse_x)

% returns positive and negative predictive values
%
% [PV,thr]=shiStatPredictiveValue(z,x)
% [PV,thr]=shiStatPredictiveValue(z,x,n_Point)
% [PV,thr]=shiStatPredictiveValue(z,x,n_Point,do_Reverse_x)
%   computes predictive values of x on z
% 
%   z             - outcome variable, >0 or true for positive outcome, <=0 
%                   or false for negative outcome
%   x             - predictor
%   n_Point       - number of thresholding points for x (default = 500), or 
%                   0 to use raw x values as thresholds
%   do_Reverse_x  - set as true if assuming negative dependence between x
%                   and z (default = false)
%   PV            - two columns, positive and negative predictive values at
%                   different threshold
%   thr           - thresholds
% 
% 
% by Zhenhao Shi @ 2015-5-18
% 

if nargin < 3
    n_Point = 500;
end

if nargin < 4
    do_Reverse_x = mean(x(z>0)) < mean(x(z<0));
end

if n_Point == 0
    thr = sort(x);
else
    thr = min(x):((max(x)-min(x))/(n_Point+2)):max(x);
    thr = thr';
end
thr = thr(2:end-2);

ppv = nan(size(thr));
npv = nan(size(thr));

if do_Reverse_x
    for i = 1:length(thr)
        pos_call = sum(x<thr(i));
        neg_call = sum(x>thr(i));
        pos_true = sum(x<thr(i) & z>0);
        neg_true = sum(x>thr(i) & z<=0);
        ppv(i) = pos_true/pos_call;
        npv(i) = neg_true/neg_call;
    end
else
    for i = 1:length(thr)
        pos_call = sum(x>thr(i));
        neg_call = sum(x<thr(i));
        pos_true = sum(x>thr(i) & z>0);
        neg_true = sum(x<thr(i) & z<=0);
        ppv(i) = pos_true/pos_call;
        npv(i) = neg_true/neg_call;
    end
end

PV = [ppv,npv];
plot(thr,ppv,'r',thr,npv,'b','LineWidth',2);
axis([thr(1) thr(end)  min(PV(:))-0.05 1.05]);
% if do_Reverse_x
%     set(gca,'XDir','reverse');
% end
legend({'PPV';'NPV'},'Location','South');
