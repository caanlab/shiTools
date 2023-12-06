function [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img,AtlasName,NodeValue,LabelListTxtFile,SummFunc)

% extracts summary value(s) of ROIs defined by an atlas (see shiMisc folder)
%
% example: [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img,'AAL',1:10)
% example: [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img,'AAL',1:10,'shiSpmTemplate_AAL_label.txt','mean')
%
% zhenhao shi



if ismember(AtlasName,shiSpmAnatLabel) && (~exist('NodeList','var') || isempty(LabelListTxtFile))
    LabelListTxtFile = which(['shiSpmTemplate_',AtlasName,'_Label.txt']);
    AtlasName0 = which(['shiSpmTemplate_',AtlasName,'.img']);
    if isempty(AtlasName0)
        AtlasName = which(['shiSpmTemplate_',AtlasName,'.nii']);
    else
        AtlasName = AtlasName0;
    end
end

%% read image

File_Img = char(Img);
for i = size(File_Img,1)
    if ~exist(deblank(File_Img(i,:)),'file')
        error('cannot find %s',File_Img(i,:))
    end
end


V_Img = spm_vol(File_Img);
[Y_Img,~] = spm_read_vols(V_Img);

n_Img = numel(V_Img);
if numel(V_Img) == 1 && length(V_Img.dim) > 3
    n_Img = V_Img.dim(4);
end


Mat_Img = V_Img(1).mat;
Dim_Img = V_Img(1).dim(1:3);

Y_Img = reshape(Y_Img,[prod(Dim_Img),n_Img])';


if ~spm_type(V_Img(1).dt(1),'nanrep')
    Y_Img(Y_Img==0) = NaN;
end



%% read node info

V_Roi = spm_vol(AtlasName);
Y_Roi = spm_read_vols(V_Roi);

if ~exist('NodeValue','var') || isempty(NodeValue)
    if exist('LabelListTxtFile','var') && ~isempty(LabelListTxtFile)
        try
            fid=fopen(LabelListTxtFile);
            NodeValue = textscan(fid,'%f%s%*[^\n]');
            fclose(fid);
            NodeValue = NodeValue{1};
        catch
            warning('cannot read node value from *Label.txt');
            NodeValue = sort(unique(Y_Roi));
            NodeValue = NodeValue(~isnan(NodeValue));
            NodeValue = NodeValue(NodeValue~=0);
        end
    end
end
NodeValue = NodeValue(:);

if exist('LabelListTxtFile','var') && ~isempty(LabelListTxtFile)
    try
        fid=fopen(LabelListTxtFile);
        LabelListTxtFile = textscan(fid,'%f%s%*[^\n]');
        fclose(fid);
        NodeName = cell(size(NodeValue));
        for i = 1:length(NodeName)
            NodeName{i} = LabelListTxtFile{2}{LabelListTxtFile{1}==NodeValue(i)};
        end
    catch
        NodeName = shiStrConcat('Region_',NodeValue);
        warning('cannot read node name from *Label.txt');
    end
end



%% compute masks

X = NaN(n_Img,length(NodeValue));

Mat_Roi = V_Roi.mat;


msk = cell(length(NodeValue),1);
for i = 1:length(NodeValue)
    msk{i} = true(Dim_Img);
end

for j = 1:Dim_Img(3)

    x_Slice_Img = spm_matrix([0 0 j]); % matrix, for a slice of Image

    Mat_Roi_to_ImgSlice = Mat_Img\Mat_Roi\x_Slice_Img; % matrix, from roi to a slice of Image
    img_slice = spm_slice_vol(V_Roi,Mat_Roi_to_ImgSlice,Dim_Img(1:2),[0 NaN]); % a slice of ROI in Image space

    for i = 1:length(NodeValue)
        msk{i}(:,:,j) = img_slice == NodeValue(i);
    end

end



%% apply masks

for i = 1:length(NodeValue)

    msk{i} = reshape(msk{i},prod(Dim_Img),1)';
    xY_Img = Y_Img(:,msk{i});


    if ~exist('SummFunc','var') || isempty(SummFunc) || strcmpi(SummFunc(1:3),'mea')
        X(:,i) = nanmean(xY_Img,2);
    elseif strcmpi(SummFunc(1:3),'med')
        X(:,i) = nanmedian(xY_Img,2);
    elseif strcmpi(SummFunc(1:3),'eig') % see spm_regions
        try
            not_nan = all(~isnan(xY_Img));
            xY_Img = xY_Img(:,not_nan);
            [m,n]   = size(xY_Img);
            if m > n
                [~,s,v] = svd(xY_Img'*xY_Img);
                s       = diag(s);
                v       = v(:,1);
                u       = xY_Img*v/sqrt(s(1));
            else
                [~,s,u] = svd(xY_Img*xY_Img');
                s       = diag(s);
                u       = u(:,1);
                v       = xY_Img'*u/sqrt(s(1));
            end
            d       = sign(sum(v));
            u       = u*d;
            X(:,i)       = u*sqrt(s(1)/n);
        catch
            X(:,i) = NaN;
        end
    else
        warning('unknown summarizing methods. use ''mean'' instead');
        X(:,i) = nanmean(xY_Img,2);
    end

end
