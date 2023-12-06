function [Estimability,EstimabilityVector] = shiSpmEstimability(Des)

% returns the estimability of an fMRI design
%
% [Estimability,EstimabilityVector] = shiSpmEstimability(Des)
% 
%   Des - Design matrix, or SPM.mat (SPM struct or file name)
%   Estimability       - Estimability of the design
%   EstimabilityVector - Estimability of each column
% 
%    ###########
% by Zhenhao Shi @ 2015-1-5
%    ###########
% 

if isstruct(Des)
    Des = Des.xX.X;
elseif exist(Des,'file') && shiStrIncl(Des,'SPM.mat$',0)
    Des = load(Des);
    Des = Des.xX.X;
end

EstimabilityVector = spm_SpUtil('isCon',Des);
Estimability = all(EstimabilityVector);

