function shiSpmR2Z(RImg,ZImg)
 
% uses Fisher-transform from correlation coefficient R image to Z image
% 
% shiSpmR2Z(RImg,ZImg)
%   RImg - string or cell array of strings for input R image(s)
%   ZImg - string or cell array of strings for output Z image(s)
% 
%    ###########
% by Zhenhao Shi @ 2015-1-6
%    ###########
% 

RImg = cellstr(char(RImg));
ZImg = cellstr(char(ZImg));

if numel(ZImg) ~= numel(RImg)
    error('unmatched file number');
end


Expression = 'atanh(i1)';

for i = 1:length(RImg)
    spm_imcalc(RImg{i},ZImg{i},Expression);
end
