function shiSpmT2R(TImg,RImg,n)

% converts student statistic T map to correlation coefficient R map
% 
% shiSpmT2R(TImg,RImg,n)
% 
%   TImg - string or cell array of strings for input T image(s)
%   RImg - string or cell array of strings for output R image(s)
%   n    - number of observations
% 
%    ###########
% by Zhenhao Shi @ 2023-5-26
%    ###########
% 

RImg = cellstr(char(RImg));
TImg = cellstr(char(TImg));

if numel(TImg) ~= numel(RImg)
    error('unmatched file number');
end

N = num2str(n);

Expression = ['sign(i1).*sqrt((i1.*i1)./(i1.*i1+(',N,'-2)))'];

for i = 1:length(TImg)
    spm_imcalc(TImg{i},RImg{i},Expression);
end
