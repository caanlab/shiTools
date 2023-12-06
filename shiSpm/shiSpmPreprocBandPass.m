function outImg = shiSpmPreprocBandPass(Img,TR,HighCutoff,LowCutoff,Prefix,existAction)

% performs band-pass filtering
%
% shiSpmPreprocBandPass(Img,TR)
% shiSpmPreprocBandPass(Img,TR,HighCutoff,LowCutoff)
% shiSpmPreprocBandPass(Img,TR,HighCutoff,LowCutoff,Prefix)
% 
%   Img         - raw images
%   TR          - repetition time (in second)
%   HighCutoff  - bandpass filter high cutoff in Hz (default 0.01)
%   LowCutoff   - bandpass filter low cutoff in Hz (default 0.08)
%
% Modified from Jin-hui Wang's gretna_preprocessing_BandPassFilter
% wjhmjy@gmail.com
%
% Zhenhao Shi, 2019-10-31
%

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'f';
end
if ~exist('HighCutoff','var') || isempty(HighCutoff)
    HighCutoff = 0.01;
end
if ~exist('LowCutoff','var') || isempty(LowCutoff)
    LowCutoff = 0.08;
end

[HighCutoff,LowCutoff] = deal(min(HighCutoff,LowCutoff),max(HighCutoff,LowCutoff));

Img = cellstr(char(Img));
[pth,nme,ext] = shiFileParts(Img);
outImg = shiStrConcat(pth,filesep,Prefix,nme,ext);

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end
if exist(outImg{1},'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',outImg{1}), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end

V = spm_vol(Img{1});
Dimension = V.dim;
Img_char = char(Img);



%Arrange the border of filtered frequency band
sampleFreq 	  = 1/TR;
sampleLength  = length(Img);
paddedLength  = shiSpmBandPass_nextpow2_one35(sampleLength);
freqPrecision = sampleFreq/paddedLength;

%-------------------------------------------------------------
%| Generate low- and high- pass filter mask (refer to
%| rest_bandpass.m written by Song Xiaowei)
%-------------------------------------------------------------
% Generate low-pass filter mask
MaskData = ones(Dimension)==1;
maskLowCutoff  = repmat(MaskData, [1, 1, 1, paddedLength]);
maskHighCutoff = maskLowCutoff;
if (LowCutoff >= (sampleFreq/2)) || (LowCutoff==0)
    maskLowCutoff(:,:,:,:) = 1;	% All pass
elseif (LowCutoff > 0) && (LowCutoff < freqPrecision)
    maskLowCutoff(:,:,:,:) = 0;	% All stop
else
    % Low pass, e.g., freq < 0.08 Hz
    idxCutoff	= round(LowCutoff * paddedLength * TR);
    idxCutoff2	= paddedLength+2 - idxCutoff;
    maskLowCutoff(:,:,:,idxCutoff+1:idxCutoff2-1) = 0;
end

% Generate high-pass filter mask
if (HighCutoff < freqPrecision)
    maskHighCutoff(:,:,:,:) = 1;	% All pass
elseif (HighCutoff >= (sampleFreq/2))
    maskHighCutoff(:,:,:,:) = 0;	% All stop
else
    % high pass, e.g., freq > 0.01 Hz
    idxCutoff	= round(HighCutoff * paddedLength *TR);
    idxCutoff2	= paddedLength+2 - idxCutoff;
    maskHighCutoff(:,:,:,1:idxCutoff-1) = 0;
    maskHighCutoff(:,:,:,idxCutoff2+1:paddedLength) = 0;
end

% filtering the images individually


FileFid = spm_vol(Img_char);
n_Img = size(Img_char,1);
nslices = FileFid(1).dim(3);
nslices_t = FileFid(1).dim(3);
if ( nslices_t ~= nslices )
    fprintf('Number of slices differ! %d %\n', nimg);
else

    Vout = FileFid;
    for k = 1:n_Img
        Vout(k).fname  = outImg{k};
        Vout(k).dt = [16 spm_platform('bigend')];
    end

    Vout = spm_create_vol(Vout);%,'noopen');
    slices = zeros([Vout(1).dim(1:2) n_Img]);
    for k = 1:nslices
        B  = spm_matrix([0 0 k]);
        for m = 1:n_Img
            slices(:,:,m) = spm_slice_vol(FileFid(m),B,FileFid(1).dim(1:2),1);
        end
        %Save the linear trend
        theTrend_Intercept = slices(:,:,1);
        theTrend_Slope = (slices(:,:,end) - theTrend_Intercept) / (sampleLength-1);
        for y = 1:sampleLength
            %remove the linear trend first
            slices(:,:,y) = slices(:,:, y) - (theTrend_Intercept + y * theTrend_Slope);
        end
        %FFT
        slicesfreq = fft(slices, paddedLength, 3);
        %Mask redundant frequency components
        FilterMask = squeeze(maskLowCutoff(:,:,k,:));
        slicesfreq(~FilterMask) = 0;
        FilterMask = squeeze(maskHighCutoff(:,:,k,:));
        slicesfreq(~FilterMask) = 0;
        %inverse FFT
        FilteredSlices = ifft(slicesfreq, paddedLength, 3);
        FilteredSlices = FilteredSlices(:,:,1:sampleLength);%remove the padded parts
        %retrend the time course

        for y = 1:n_Img
            %add the linear trend after filter
            FilteredSlices(:,:, y) = FilteredSlices(:,:, y) + (theTrend_Intercept+y*theTrend_Slope);
        end

        % write out the slice for all volumes
        for p = 1:n_Img
            Vout(p) = spm_write_plane(Vout(p),FilteredSlices(:,:,p),k);
        end
    end
    %         Vout = spm_close_vol(Vout);
end


function Result = shiSpmBandPass_nextpow2_one35(n)
%Compute the min length for FFT according to AFNI's algorithm, By Xiao-Wei Song
%------------------------------------------------------------------------------------------------------------------------------
%	Copyright(c) 2007~2010
%	State Key Laboratory of Cognitive Neuroscience and Learning in Beijing Normal University
%	Written by Xiao-Wei Song 
%	http://resting-fmri.sourceforge.net
% 	<a href="Dawnwei.Song@gmail.com">Mail to Author</a>: Xiaowei Song
%	Version=1.0;
%	Release=20070903;

    if length(n) > 1
        n = cast(length(n),class(n));
    end
    if n < 16
        Result = 2^nextpow2(n);
        return;
    end 
    
    limit = nextpow2(n);             % n = 134, limit = 8
    tbl = 2^(limit-1):2^limit;       % tbl = 128, 129, ... , 256
    tbl = tbl(tbl>=n);               % tbl = 134, 135, ... , 256
    for x = 1:length(tbl)
        Result = tbl(x);
        [f,~] = log2(Result);
        if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
            return;
        end
        if mod(Result,3*5) == 0        
            y = Result / (3*5);
            [f,~] = log2(y);
            if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,3) == 0        
            y = Result / 3;
            [f,~] = log2(y);
            if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,5) == 0        
            y = Result / 5;
            [f,~] = log2(y);
            if ~isempty(f) && f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
    end
    Result = NaN;    % Should not reach, except when n=1

