function shiSpmT2Z(TImg,ZImg,df)

% converts student statistic T map to unit normal Z map
% 
% shiSpmT2Z(TImg,RImg,n)
% 
%   TImg - string or cell array of strings for input T image(s)
%   ZImg - string or cell array of strings for output R image(s)
%   df   - degree of freedom
% 
%    ###########
% by Zhenhao Shi @ 2016-3-24
%    ###########
% 

ZImg = cellstr(char(ZImg));
TImg = cellstr(char(TImg));

if numel(TImg) ~= numel(ZImg)
    error('unmatched file number');
end

DF = num2str(df);

Expression = ['spm_t2z(i1,75',DF,')'];

for i = 1:length(TImg)
    spm_imcalc(TImg{i},ZImg{i},Expression);
end;
