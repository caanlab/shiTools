function shiSpmR2T(RImg,TImg,n)
 
% converts correlation coefficient R map to student statistic T map
% 
% shiSpmR2T(RImg,TImg,n)
%   RImg - string or cell array of strings for input R image(s)
%   TImg - string or cell array of strings for output T image(s)
%   n    - number of observations
% 
%    ###########
% by Zhenhao Shi @ 2018-4-21
%    ###########
% 

RImg = cellstr(char(RImg));
TImg = cellstr(char(TImg));

if numel(TImg) ~= numel(RImg)
    error('unmatched file number');
end

for i = 1:length(TImg)
    shispmr2t(RImg{i},TImg{i},n);
end


function shispmr2t(r,t,n)

rV = spm_vol(r);
rY = spm_read_vols(rV);

tV = struct('fname',   t,...
            'dim',     rV.dim(1:3),...
            'dt',      [64 spm_platform('bigend')],...
            'mat',     rV.mat,...
            'descrip', ['shiSpmR2T:SPM{T_[',num2str(n-2),']}']);
tY = rY.*sqrt((n-2)./(1-rY.^2));

spm_write_vol(tV,tY);