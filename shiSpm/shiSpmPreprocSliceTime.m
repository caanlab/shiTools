function [outImg,matlabbatch] = shiSpmPreprocSliceTime(Img,Tr,Ta,SliceOrder,RefSlice,Prefix,existAction)

% performs preprocessing: slice time correction
%
% outImg = shiSpmPreprocSliceTime(Img,Tr,Ta,SliceOrder)
% outImg = shiSpmPreprocSliceTime(Img,Tr,Ta,SliceOrder,RefSlice)
% outImg = shiSpmPreprocSliceTime(Img,Tr,Ta,SliceOrder,RefSlice,Prefix)
%
%   Img                 - cellstr, functional images
%   Tr                  - number, repetition time (in second)
%   Ta                  - number, time of acquisition (in second) (according to SPM, Ta will be unused if slice order is provided as acquisition time in millisecond)
%   SliceOrder          - slice order or slice time in millisecond
%   RefSlice            - reference slice for slice timing (spatial slice index or slice time in millisecond)
%
% Zhenhao Shi, 2019-10-30
%


Img = cellstr(char(Img));
if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'a';
end

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

if (~exist('RefSlice','var')||isempty(RefSlice)) && (~exist('SliceOrder','var')||isempty(SliceOrder))
    fprintf('\nWARNING:\n    Slice Timing: no slice info provided, skipping and simply copying images...\n\n');
    for k = 1:length(Img)
        copyfile(Img{k},outImg{k});
    end
    matlabbatch = [];
    return;
end

if ~isequal(1:numel(SliceOrder),sort(SliceOrder))
    unit = 'slice times';
else
    unit = 'slice indices';
end

nSlices = numel(SliceOrder);

if ~exist('RefSlice','var') || isempty(RefSlice) || ~(RefSlice>=0)
    switch unit
        case 'slice indices'
            RefSlice = SliceOrder(floor(median(1:nSlices)));
        case 'slice times'
            RefSlice = SliceOrder(abs(SliceOrder-median(SliceOrder))==min(abs(SliceOrder-median(SliceOrder))));
            RefSlice = RefSlice(end);
    end
end

% switch unit
%     case 'slice times'
%         timing = [0,Tr];
%     case 'slice indices'
%         timing = [Ta/(numel(SliceOrder)-1),Tr-Ta];
% end
% 
% spm_slice_timing(char(Img), SliceOrder, RefSlice, timing, Prefix);


%-----------------------------------------------------------------------
% Job saved on 30-Oct-2019 13:33:16 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.temporal.st.scans = {Img}';
%%
matlabbatch{1}.spm.temporal.st.nslices = nSlices;
matlabbatch{1}.spm.temporal.st.tr = Tr;
matlabbatch{1}.spm.temporal.st.ta = Ta;
matlabbatch{1}.spm.temporal.st.so = SliceOrder;
matlabbatch{1}.spm.temporal.st.refslice = RefSlice;
matlabbatch{1}.spm.temporal.st.prefix = Prefix;

spm_jobman('serial',matlabbatch);