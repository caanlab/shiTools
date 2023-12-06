function [X_filtered,ActualHighpass,ActualLowpass] = shiStatBandPass(X,SampRate,HighCutoff,LowCutoff)
%
% performs band-pass filtering of time series
%
% X             - vector or matrix time series data
% SampRate      - sampling rate in Hz (NOT in sec) (default = 1 Hz)
% HighCutoff    - high-pass cutoff in Hz (default = 0.01 Hz)
% LowCutoff     - low-pass cutoff in Hz (default = 0.1 Hz)
%
% Modified from Jin-hui Wang's gretna_preprocessing_BandPassFilter
% wjhmjy@gmail.com

if ~exist('HighCutoff','var') || isempty(HighCutoff)
    HighCutoff = 0.01;
end
if ~exist('LowCutoff','var') || isempty(LowCutoff)
    LowCutoff = 0.1;
end
if ~exist('SampRate','var') || isempty(SampRate)
    SampRate = 1;
end

[HighCutoff,LowCutoff] = deal(min(HighCutoff,LowCutoff),max(HighCutoff,LowCutoff));

if size(X,2)>1
    X_filtered = nan(size(X));
    [ActualHighpass,ActualLowpass] = deal(nan(1,size(X,2)));
    for i = 1:size(X,2)
        [X_filtered(:,i),ActualHighpass(:,i),ActualLowpass(:,i)] = shiBandPass(X(:,i),SampRate,HighCutoff,LowCutoff);
    end
    return;
end

sampleLength = length(X);
X = reshape(X,length(X),1);
paddedLength = filter_nextpow2_one35(sampleLength);
freqPrecision= SampRate/paddedLength;

%%

maskLowPass =ones(paddedLength,1);
maskHighPass=maskLowPass;

if (LowCutoff>=(SampRate/2))||(LowCutoff==0)
    maskLowPass(:)=1;	%All pass
    ActualLowpass = SampRate/2;
elseif (LowCutoff>0)&&(LowCutoff< freqPrecision)
    maskLowPass(:)=0;	% All stop
    ActualLowpass = 0;
else
    % Low pass, e.g., freq < 0.08 Hz
    idxCutoff	=round(LowCutoff *paddedLength /SampRate);
    ActualLowpass = idxCutoff./(1/SampRate*paddedLength);
    idxCutoff2	=paddedLength+2 -idxCutoff;
    maskLowPass(idxCutoff+1:idxCutoff2-1)=0;
end

if (HighCutoff < freqPrecision)
    maskHighPass(:)=1;	%All pass
    ActualHighpass=0;
elseif (HighCutoff >= (SampRate/2))
    maskHighPass(:)=0;	% All stop
    ActualHighpass = SampRate/2;
else
    % high pass, e.g., freq > 0.01 Hz
    idxCutoff	=round(HighCutoff *paddedLength /SampRate);
    ActualHighpass = idxCutoff./(paddedLength/SampRate);
    idxCutoff2	=paddedLength+2 -idxCutoff;
    maskHighPass(1:idxCutoff-1)=0;
    maskHighPass(idxCutoff2+1:paddedLength)=0;
end


%%

theTrend_Intercept=X(1);
theTrend_Slope= (X(end) -theTrend_Intercept)/(sampleLength-1);
for y=1:sampleLength
    X(y)=X(y) -(theTrend_Intercept + y*theTrend_Slope);
end

Xfreq =fft(X, paddedLength);

Xfreq(~maskLowPass)=0;
Xfreq(~maskHighPass)=0;

X_filtered =ifft(Xfreq, paddedLength);
X_filtered =X_filtered(1:sampleLength);

for y=1:sampleLength
    X_filtered(y)=X_filtered(y) + (theTrend_Intercept+y*theTrend_Slope);
end


%%

function Result = filter_nextpow2_one35(n)

    if length(n)>1
        n = cast(length(n),class(n));
    end
    if n<16
        Result =2^nextpow2(n);
        return;
    end 
    
    limit =nextpow2(n);             %n=134, limit=8
    tbl=2^(limit-1):2^limit;      %tbl =128, 129, ... , 256
    tbl =tbl(tbl>=n);          %tbl =134, 135, ... , 256
    for x=1:length(tbl)
        Result =tbl(x);
        [f,~]=log2(Result);
        if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
            return;
        end
        if mod(Result,3*5)==0        
            y= Result /(3*5);
            [f,~]=log2(y);
            if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,3)==0        
            y= Result /3;
            [f,~]=log2(y);
            if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,5)==0        
            y= Result /5;
            [f,~]=log2(y);
            if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
    end
    Result =NaN;    % Should not reach, except when n=1

% csfft_nextup35 in AFNI list 1~1024, 20070516, dawnsong
% 2
% 4
% 6
% 8
% 10
% 12
% 16
% 20
% 24
% 30
% 32
% 40
% 48
% 60
% 64
% 80
% 96
% 120
% 128
% 160
% 192
% 240
% 256
% 320
% 384
% 480
% 512
% 640
% 768
% 960
% 1024