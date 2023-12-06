function Yf = shiSpmFilter(Y,HighPassPeriod,TR,WindowLength)

% performs SPM-like high-pass filter on 1-D time series or .img files (calling of "shiSpmImgCalc2" has been outdated; needs further update)
%
% Yf = shiSpmFilter(Y,HighPassPeriod,TR)
% Yf = shiSpmFilter(Y,HighPassPeriod,TR,WindowLength)
% 
%   Y               - raw 1-D time series or .img file names
%   HighPassPeriod  - high-pass cutoff given by period in second (e.g. 128)
%   TR              - sampling rate in second (e.g. 2)
%   WindowLength    - length of each sub-process if any (default = the
%                     full length of Y). E.g., if there are two runs, each
%                     contains 150 TRs, then the length of Y must be 300,
%                     and WindowLength should be given as [150,150].
%   Yf              - filtered time series or 'filtered_*.img' files
% 
% requires
%   SPM 8
%   shiSpmImgCalc
% 
%    ###########
% by Zhenhao Shi @ 2015-2-12
%    ###########
% 

if isnumeric(Y)
    if length(size(Y))~=2 || min(size(Y)) ~=1
        error('Y must be a 1-D vector');
    end
else
    Y = cellstr(char(Y));
    for i = 1:length(Y)
        yyy = shiFullFileName(Y{i});
        Y{i} = yyy{1};
        [xpath,xname,xext] = fileparts(Y{i});
        Yf{i,1} = fullfile(xpath,['filtered_',xname,xext]);
    end;
    Script = ['Output = shiSpmFilter(Input,',num2str(HighPassPeriod),',',num2str(TR),');'];
    shiSpmImgCalc2(Y,Script,Yf);
    return;
end

Y = Y(:);

if nargin < 4
    WindowLength = length(Y);
end

WindowLength = WindowLength(:);
if ~isequal(sum(WindowLength),length(Y))
    error('sum of window length must equal to length of Y');
end

WindowLength = WindowLength(WindowLength~=0);

for i = 1:numel(WindowLength)
    K(i).RT = TR;
    K(i).row = sum(WindowLength(1:i-1))+1:sum(WindowLength(1:i));
    K(i).HParam = HighPassPeriod;
end;

K = spm_filter(K);

Yf = spm_filter(K,Y);