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
% by Zhenhao Shi @ 2025.1.19
%    ###########
% 


TImg = cellstr(char(TImg));
ZImg = cellstr(char(ZImg));

if numel(ZImg) ~= numel(TImg)
    error('unmatched file number');
end

for i = 1:length(ZImg)
    shispmt2z(TImg{i},ZImg{i},df);
end


function shispmt2z(t,z,df)

tV = spm_vol(t);
tY = spm_read_vols(tV);

zV = struct('fname',   z,...
            'dim',     tV.dim(1:3),...
            'dt',      [64 spm_platform('bigend')],...
            'mat',     tV.mat,...
            'descrip', 'shiSpmT2Z');
zY = spm_t2z(tY,df,75);

spm_write_vol(zV,zY);