function [Roi,RoiFormat,RoiName] = shiSpmRoiFormat(Roi,Radius)

% formats ROI into a readable cell array
%
% Input: Roi can be:
%   - a nx3 matrix for spherical ROIs, each row being MNI x/y/z, default radius = 5 mm
%   - a nx4 matrix for spherical ROIs, each row being MNI x/y/z plus radius (mm)
%   - a char matrix for binary masks, each row being a mask filename
%   - a single struct for ROIs defined by one or more atlas, with the following fields:
%       - .atlas, one or n strings, each being an atlas filename or the name of a built-in atlas (e.g. 'AAL', 'HOA'. see "shiMisc" folder)
%       - .value, scalar or nx1 vector indicating the values corresponding to the parcels of interest; may also be string(s) of parcel names, if .list exist (see below)
%       - .list (optional), one or n strings, each being a .txt file that contains a column of parcel valus and another column of parcel names
%       (all fields must have the same length of n, unless when .atlas or .list have length of 1, in which case they will be applied to all values in .value)
%   - a struct array similar to the above, where each struct has one string in .atlas, one value (numeric or string) in .value, and optionally one filename in .list
%   - a cell array, each cell defining either (i) one spherical ROI, (ii) one binary mask, (iii) or one parcel in one atlas
% Output: Roi will be arranged as a cell array described above; RoiFormat = 'Atlas', 'Sphere', or 'Mask'
%
% zhenhao shi 2018-6-19
%


if ~exist('Radius','var') || isempty(Radius)
    Radius = 5; % not used unless not specified otherwise
end

%% convert vector, string, and struct to cell

Roi_orig = Roi;
IS_NUM = false;
IS_CHAR = false;
IS_STRUCT = false;
if isnumeric(Roi_orig) % n-by-4 coordinate-radius matrix, or n-by-3 coordinate matrix
    IS_NUM = true;
    if size(Roi_orig,2) == 3
        Roi_orig(:,end+1) = Radius;
    elseif size(Roi_orig,2) ~= 3 && size(Roi_orig,2) ~= 4
        error('Roi in matrix format: must be n-by-4 (x,y,z,radius) or n-by-3 (x,y,z) matrix');
    end
    Roi = mat2cell(Roi_orig,ones(size(Roi_orig,1),1),size(Roi_orig,2));
elseif ischar(Roi_orig) % character matrix for single-label ROI image(s)
    IS_CHAR = true;
    Roi = cellstr(Roi_orig);
    for i = 1:length(Roi)
        if ~exist(Roi{i},'file')
            error('Roi # %d: cannot find file %s',i,Roi{i});
        end
    end
elseif isstruct(Roi_orig)
    IS_STRUCT = true;
    if ~isfield(Roi_orig,'atlas') || ~isfield(Roi_orig,'value')
        error('Roi in struct format: must have fields .atlas and .value');
    end
    if numel(Roi_orig) > 1 % struct array with a string field .atlas and a scalar/string field .value
        Roi = cell(numel(Roi_orig),1);
        if ~isfield(Roi_orig,'list')
            [Roi_orig.list] = deal('');
        end
        for i = 1:numel(Roi_orig)
            [Roi{i}.atlas,Roi{i}.value,Roi{i}.list,Roi{i}.label] = shiSpmRoiFormat_parcel(i,Roi_orig(i).atlas,Roi_orig(i).value,Roi_orig(i).list);
        end
    else % single struct with a string/cellstr field .atlas and a vector/cellstr field .value
        if ischar(Roi_orig.value)
            Roi_orig.value = cellstr(Roi_orig.value);
        end
        if ~isfield(Roi_orig,'list')
            Roi_orig.list = repmat({''},numel(Roi_orig.value),1);
        else
            Roi_orig.list = cellstr(char(Roi_orig.list));
            if length(Roi_orig.list) == 1
                Roi_orig.list = repmat(Roi_orig.list,numel(Roi_orig.value),1);
            end
        end
        Roi_orig.atlas = cellstr(char(Roi_orig.atlas));
        if length(Roi_orig.atlas) == 1
            Roi_orig.atlas = repmat(Roi_orig.atlas,numel(Roi_orig.value),1);
        end
        if length(Roi_orig.list) ~= length(Roi_orig.atlas) || length(Roi_orig.atlas) ~= length(Roi_orig.value)
            error('Roi in struct format: make sure that .atlas, .value, and, if existing, .list, all have the same length');
        end
        Roi = cell(numel(Roi_orig.value),1);
        for i = 1:numel(Roi_orig.value)
            [Roi{i}.atlas,Roi{i}.value,Roi{i}.list,Roi{i}.label] = shiSpmRoiFormat_parcel(i,Roi_orig.atlas(i),Roi_orig.value(i),Roi_orig.list(i));
        end
    end
elseif ~iscell(Roi_orig)
    error('Roi: unknown ROI format');
end


%% check each cell of ROI

if IS_STRUCT
    RoiFormat = repmat({'Atlas'},length(Roi),1);
elseif IS_NUM
    RoiFormat = repmat({'Sphere'},length(Roi),1);
elseif IS_CHAR
    RoiFormat = repmat({'Mask'},length(Roi),1);
else
    RoiFormat = cell(length(Roi),1);
    for r = 1:length(Roi)
        if isnumeric(Roi{r})
            if ~isequal(size(Roi{r}),[1,3]) && ~isequal(size(Roi{r}),[1,4])
                error('Roi # %d: coordinate information must be [x,y,z,radius] or [x,y,z] (default radius = 5)',r);
            elseif size(Roi{r},2)==3
                Roi{r} = [Roi{r},Radius];
            elseif ~(Roi{r}(4)>=0)
                error('Roi # %d: radius must be >= 0 (when =0, extract nearest voxel)',r);
            end
            RoiFormat{r} = 'Sphere';
        elseif ischar(Roi{r})
            if ~exist(Roi{r},'file')
                error('Roi # %d: cannot find file %s',r,Roi{r});
            end
            RoiFormat{r} = 'Mask';
        elseif isstruct(Roi{r}) && all(isfield(Roi{r},{'atlas','value'}))
            RoiFormat{r} = 'Atlas';
            if ~isfield(Roi{r},'list')
                Roi{r}.list = [];
            end
            [Roi{r}.atlas,Roi{r}.value,Roi{r}.list,Roi{r}.label] = shiSpmRoiFormat_parcel(r,Roi{r}.atlas,Roi{r}.value,Roi{r}.list);
        elseif isstruct(Roi{r})
            error('Roi # %d: must have fields .atlas and .value',r);
        else
            error('Roi # %d: unknown ROI format',r);
        end
    end
end
Roi = Roi(:);

RoiName = cell(size(Roi));
for r = 1:length(Roi)
    switch RoiFormat{r}
        case 'Sphere'
            RoiName{r} = sprintf('%.4g_%.4g_%.4g_%.4g',Roi{r});
        case 'Mask'
            [~,RoiName{r}] = fileparts(Roi{r});
        case 'Atlas'
            if ~isfield(Roi{r},'label') || ~ischar(Roi{r}.label) || ~isempty(Roi{r}.label)
                RoiName{r} = [Roi{r}.atlas,'__',sprintf('%.5g',Roi{r}.value)];
            else
                RoiName{r} = [Roi{r}.atlas,'__',Roi{r}.label];
            end
    end
end

function [atlas,value,list,label] = shiSpmRoiFormat_parcel(ind,atlas,value,list)

label = '';
atlas = char(atlas);
if exist(atlas,'file')
    atlas_type = 'file';
elseif ismember(atlas,shiSpmAnatLabel)
    atlas_type = 'preset';
else
    error('Roi # %d: cannot find atlas %s',ind,atlas);
end

if iscellstr(value) || ischar(value)
    value = char(value);
    value_type = 'char';
    if size(value,1) ~= 1
        error('Roi # %d: can only process one parcel at a time',ind);
    end
elseif isnumeric(value)
    value_type = 'num';
    if length(value) ~= 1
        error('Roi # %d: can only process one parcel at a time',ind);
    end
else
    error('Roi # %d: parcel value can only be numeric or string',ind);
end

list = char(list);
if exist(list,'file')
    list_type = 'txt';
else
    list_type = 'none';
end

switch [atlas_type,'_',list_type]
    case 'file_none'
        if isequal(value_type,'char')
            error('Roi # %d: no label .txt file found, .parcel has to be a number',ind);
        elseif isequal(value_type,'num')
            label = '';
        end
    case 'preset_txt'
        atlas0 = atlas;
        atlas = which(['shiSpmTemplate_',atlas0,'.nii']);
        if isempty(atlas)
            atlas = which(['shiSpmTemplate_',atlas0,'.img']);
        end
    case 'preset_none'
        atlas0 = atlas;
        list = which(['shiSpmTemplate_',atlas0,'_Label.txt']);
        atlas = which(['shiSpmTemplate_',atlas0,'.nii']);
        if isempty(atlas)
            atlas = which(['shiSpmTemplate_',atlas0,'.img']);
        end
        list_type = 'txt';
end

if isequal(value_type,'num')
    if isequal(list_type,'txt')
        fid = fopen(list);
        LabelList=textscan(fid,'%d%s%*[^\n]');
        fclose(fid);
        LabelNum = LabelList{1};
        LabelList = LabelList{2};
        label = LabelList{LabelNum==value};
        if isempty(label)
            error('Roi # %d: cannot find parcel value=%d',ind,value);
        end
    end
elseif isequal(value_type,'char')
    label = value;
    fid = fopen(list);
    LabelList=textscan(fid,'%d%s%*[^\n]');
    fclose(fid);
    LabelNum = LabelList{1};
    LabelList = LabelList{2};
    value = LabelNum(strcmp(label,LabelList));
    if isempty(value)
        error('Roi # %d: cannot find parcel value=%d',ind,value);
    end
end


