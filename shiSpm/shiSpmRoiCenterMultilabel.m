function [XYZ,NodeValue,NodeName] = shiSpmRoiCenterMultilabel(AtlasName,NodeValue,LabelListTxtFile)

% returns the centers of mass of parcels
%
% example: [X,NodeValue,NodeName] = shiSpmRoiCenterMultilabel('AAL')
% example: [X,NodeValue,NodeName] = shiSpmRoiCenterMultilabel('AAL',1:10)
% example: [X,NodeValue,NodeName] = shiSpmRoiCenterMultilabel('AAL',1:10,'shiSpmTemplate_AAL_label.txt')
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



%% read node info

V = spm_vol(AtlasName);
[Y,xyz] = spm_read_vols(V);

if ~exist('NodeValue','var') || isempty(NodeValue)
    if exist('LabelListTxtFile','var') && ~isempty(LabelListTxtFile)
        try
            fid=fopen(LabelListTxtFile);
            NodeValue = textscan(fid,'%f%s%*[^\n]');
            fclose(fid);
            NodeValue = NodeValue{1};
        catch
            warning('cannot read node value from *Label.txt');
            NodeValue = sort(unique(Y));
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
        if isnumeric(NodeValue)
            NodeName = cell(size(NodeValue));
            for i = 1:length(NodeName)
                NodeName{i} = LabelListTxtFile{2}{LabelListTxtFile{1}==NodeValue(i)};
            end
        else
            NodeName = cellstr(NodeValue);
            NodeValue = nan(size(NodeName));
            for i = 1:numel(NodeName)
                NodeValue(i) = LabelListTxtFile{1}(matches(LabelListTxtFile{2},NodeName{i}));
            end
        end
    catch
        NodeName = shiStrConcat('Region_',NodeValue);
        warning('cannot read node name from *Label.txt');
    end
end

%%

XYZ = nan(length(NodeValue),3);
for i = 1:length(NodeValue)
    XYZ(i,:) = mean(xyz(:,Y(:)==NodeValue(i)),2)';
end