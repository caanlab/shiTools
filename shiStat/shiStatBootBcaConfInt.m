function CI = shiStatBootBcaConfInt(theta, theta_boot, theta_jack, alpha)

% computes BCa confidence interval for a parameter
%
% Syntax:
%   CI = shiRBootBcaConfInt(theta, theta_boot, theta_jack, alpha)
%
% Input:
%   theta: scalar value of parameter estimated from original data
%   theta_boot: vector of parameter estimates from bootstrap
%   theta_jack: vector of parameter estimates from jackknife
%   alpha: alpha level (default=0.05, i.e. 95% BCa CI)
%
% Output:
%   CI: BCa confidence interval
%
% Reference:
%   Efron, B. (1987). Better bootstrap confidence intervals. Journal of the American statistical Association, 82(397), 171-185.

if nargin < 4
    alpha = 0.05; % default alpha level
end

% Calculate the jackknife difference
theta_jack_diff = mean(theta_jack) - theta_jack;

% Calculate acceleration value (a)
a = sum(theta_jack_diff.^3) / (6 * sum(theta_jack_diff.^2) ^ 1.5);

% Calculate z0
z0 = norminv(mean(theta_boot <= theta));

% Calculate z_alpha1 and z_alpha2
z_alpha1 = norminv(alpha/2);
z_alpha2 = norminv(1-alpha/2);

% Calculate alpha1 and alpha2
alpha1 = normcdf(z0 + (z0 + z_alpha1) / (1 - a * (z0 + z_alpha1)));
alpha2 = normcdf(z0 + (z0 + z_alpha2) / (1 - a * (z0 + z_alpha2)));

% Compute BCa confidence interval
CI = quantile(theta_boot, [alpha1, alpha2]);
