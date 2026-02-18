function [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img,AtlasName,NodeValue,LabelListTxtFile,SummFunc)

% extracts summary value(s) of ROIs defined by an atlas (see shiSpmAnatLabel)
%
% example: [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img,'AAL',1:10)
% example: [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img,{'AAL','HOA'},1:10) (will return output in cells)
% example: [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img,'AAL',1:10,'shiSpmTemplate_AAL_label.txt','mean')
%
% zhenhao shi 2026.2.18



%% determine cell vs. matrix output for X

if iscell(AtlasName) || ( ischar(AtlasName) && ~isscalar(cellstr(AtlasName)) )
    IS_CELL = true;
else
    IS_CELL = false;
end
AtlasName = cellstr(AtlasName);
if ~IS_CELL, assert(isscalar(AtlasName)); end

[X,NodeName] = deal(cell(size(AtlasName)));


%% check input

File_Img = char(Img);
for k = size(File_Img,1)
    if ~exist(deblank(File_Img(k,:)),'file')
        error('cannot find %s',File_Img(k,:));
    end
end

if ~exist('NodeValue','var') || isempty(NodeValue)
    NodeValue = repmat({''},size(AtlasName));
else
    NodeValue = cellstr(NodeValue);
    if isscalar(NodeValue)
        NodeValue = repmat(NodeValue,size(AtlasName));
    end
end
assert(numel(NodeValue) == numel(AtlasName));

if ~exist('LabelListTxtFile','var') || isempty(LabelListTxtFile)
    LabelListTxtFile = repmat({''},size(AtlasName));
else
    LabelListTxtFile = cellstr(LabelListTxtFile);
    if isscalar(LabelListTxtFile)
        LabelListTxtFile = repmat(LabelListTxtFile,size(AtlasName));
    end
end
assert(numel(LabelListTxtFile) == numel(AtlasName));

if ~exist('SummFunc','var') || isempty(SummFunc), SummFunc = 'mea'; end


%% chunk

ChkSz = 500;
if size(File_Img,1) > ChkSz && ~strcmpi(SummFunc(1:3),'eig') % no chunking if eig
    s = sprintf('RoiXtr chunk: 0/%d',ceil(size(File_Img,1)./ChkSz));
    fprintf('%s',s);
    [xX,xNodeValue,xNodeName] = deal(cell([ceil(size(File_Img,1)./ChkSz),1]));
    for i = 1:ceil(size(File_Img,1)./ChkSz)
        fprintf(repmat('\b',1,length(s)));
        s = sprintf('RoiXtr chunk: %d/%d',i,ceil(size(File_Img,1)./ChkSz));
        fprintf('%s',s);
        [xX{i},xNodeValue{i},xNodeName{i}] = shiSpmRoiXtrMultilabel( Img( (i*ChkSz-ChkSz+1):min(i*ChkSz,size(File_Img,1)) ,:), AtlasName,NodeValue,LabelListTxtFile,SummFunc);
    end
    for a = 1:numel(AtlasName)
        for i = 1:length(xX)
            X{a} = [X{a};xX{i}{a}];
        end
    end
    NodeValue = xNodeValue{1};
    NodeName = xNodeName{1};
    if ~IS_CELL
        X = X{1};
        NodeValue = NodeValue{1};
        NodeName = NodeName{1};
    end

    fprintf('\n');
    return;
end


%% determine summary functions

if isempty(SummFunc) || strcmpi(SummFunc(1:3),'mea')
    SUMMFUNC = @shi_roixtr_mean;
elseif strcmpi(SummFunc(1:3),'med')
    SUMMFUNC = @shi_roixtr_median;
elseif strcmpi(SummFunc(1:3),'eig') % see spm_regions
    SUMMFUNC = @shi_roixtr_eig;
else
    warning('unknown summarizing methods. use ''mean'' instead');
    SUMMFUNC = @shi_roixtr_mean;
end


%% extract images

V_Img = spm_vol(File_Img);
[Y_Img,~] = spm_read_vols(V_Img);

n_Img = numel(V_Img);
if isscalar(V_Img) && length(V_Img.dim) > 3
    n_Img = V_Img.dim(4);
end

Mat_Img = V_Img(1).mat;
Dim_Img = V_Img(1).dim(1:3);

Y_Img = reshape(Y_Img,[prod(Dim_Img),n_Img])';

if ~spm_type(V_Img(1).dt(1),'nanrep')
    Y_Img(Y_Img==0) = NaN;
end


%% main loop

for i = 1:numel(AtlasName)

    %% match existing atlas

    if ismember(AtlasName{i},shiSpmAnatLabel) && isempty(LabelListTxtFile{i})
        LabelListTxtFile{i} = which(['shiSpmTemplate_',AtlasName{i},'_Label.txt']);
        AtlasName{i} = which(['shiSpmTemplate_',AtlasName{i},'.nii']);
    end

    %% read node info

    V_Roi = spm_vol(AtlasName{i});
    Y_Roi = spm_read_vols(V_Roi);

    if isempty(NodeValue{i})
        if ~isempty(LabelListTxtFile{i})
            try
                fid=fopen(LabelListTxtFile{i});
                NodeValue{i} = textscan(fid,'%f%s%*[^\n]');
                fclose(fid);
                NodeValue{i} = NodeValue{i}{1};
            catch
                warning('cannot read node value from *Label.txt');
                NodeValue{i} = sort(unique(Y_Roi));
                NodeValue{i} = NodeValue{i}(~isnan(NodeValue{i}));
                NodeValue{i} = NodeValue{i}(NodeValue{i}~=0);
            end
        end
    end
    NodeValue{i} = NodeValue{i}(:);

    if ~isempty(LabelListTxtFile{i})
        try
            fid=fopen(LabelListTxtFile{i});
            LabelListTxtFile{i} = textscan(fid,'%f%s%*[^\n]');
            fclose(fid);
            NodeName{i} = cell(size(NodeValue{i}));
            for k = 1:length(NodeName{i})
                NodeName{i}{k} = LabelListTxtFile{i}{2}{LabelListTxtFile{i}{1}==NodeValue{i}(k)};
            end
        catch
            NodeName{i} = shiStrConcat('Region_',NodeValue{i});
            warning('cannot read node name from *Label.txt');
        end
    end

    %% compute masks

    X{i} = NaN(n_Img,length(NodeValue{i}));

    Mat_Roi = V_Roi.mat;

    msk = cell(length(NodeValue{i}),1);
    for k = 1:length(NodeValue{i})
        msk{k} = true(Dim_Img);
    end

    for j = 1:Dim_Img(3)

        x_Slice_Img = spm_matrix([0 0 j]); % matrix, for a slice of Image

        Mat_Roi_to_ImgSlice = Mat_Img\Mat_Roi\x_Slice_Img; % matrix, from roi to a slice of Image
        img_slice = spm_slice_vol(V_Roi,Mat_Roi_to_ImgSlice,Dim_Img(1:2),[0 NaN]); % a slice of ROI in Image space

        for k = 1:length(NodeValue{i})
            msk{k}(:,:,j) = img_slice == NodeValue{i}(k);
        end

    end

    %% apply masks

    for k = 1:length(NodeValue{i})

        msk{k} = reshape(msk{k},prod(Dim_Img),1)';
        xY_Img = Y_Img(:,msk{k});
        X{i}(:,k) = SUMMFUNC(xY_Img);

    end

end

if ~IS_CELL
    X = X{1};
    NodeName = NodeName{1};
    NodeValue = NodeValue{1};
end


%%

function x = shi_roixtr_mean(y)
x = nanmean(y,2);

function x = shi_roixtr_median(y)
x = nanmedian(y,2);

function x = shi_roixtr_eig(y)
try
    not_nan = all(~isnan(y));
    y = y(:,not_nan);
    [m,n]   = size(y);
    if m > n
        [~,s,v] = svd(y'*y);
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
    else
        [~,s,u] = svd(y*y');
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
    end
    d       = sign(sum(v));
    u       = u*d;
    x       = u*sqrt(s(1)/n);
catch
    x = NaN;
end

