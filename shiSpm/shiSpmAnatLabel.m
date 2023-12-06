function [FirstLabel,AllLabel] = shiSpmAnatLabel(AtlasName,Roi,Radius,SupressPrint)
 
% returns the anatomical label of ROIs
%
% FirstLabel = shiSpmAnatLabel(   []    ,Roi)
% FirstLabel = shiSpmAnatLabel(AtlasName,Roi)
% FirstLabel = shiSpmAnatLabel(AtlasName,Roi,Radius)
% FirstLabel = shiSpmAnatLabel(AtlasName,Roi,Radius,SupressPrint)
% [FirstLabel,AllLabel] = shiSpmAnatLabel(AtlasName,Roi,Radius,SupressPrint)
% [Atlas_nii,AtlasLabel_Txt] = shiSpmAnatLabel(AtlasName_inShiMisc)
% 
%   AtlasName - string, 'AAL', 'AALFullName', etc
%   Roi       - see shiSpmRoiXtr
%   Radius    - a positive number in millimeter, to specify the radius of
%               spheric ROIs, but only when variable Roi contains 
%               coordinates (default = 5)
%   FirstLabel- Each row corresponds to an ROI, most likely label
%   AllLabel  - Each row corresponds to an ROI. 1st columnis percentage of
%               voxels, 2nd column is AAL labels, 3rd column is number of
%               voxels in the atlas space
% 
%    ###########
% by Zhenhao Shi @ 2018-6-14
%    ###########
% 

cwd=pwd;

cd(fullfile(shiTools,'shiMisc'));
% [~,AtlasList1] = shiFileParts(shiFileName('shiSpmTemplate_*.img'));
% [~,AtlasList2] = shiFileParts(shiFileName('shiSpmTemplate_*.nii'));
[~,AtlasList] = shiFileParts(shiFileName('shiSpmTemplate_*.nii'));
% AtlasList = [AtlasList1;AtlasList2];
AtlasList = char(AtlasList);
AtlasList = AtlasList(:,16:end);
AtlasList = cellstr(AtlasList);

cd(cwd);
% 
% AtlasList = {
%     'AAL'
%     'AALFullName'
%     'Brodmann'
%     'Dosenbach'
%     'HOA'
%     'Power264'
%     'Yeo7NetLiberal'
%     'Yeo7NetTight'
%     'Yeo17NetLiberal'
%     'Yeo17NetTight'
%     };
if nargin == 0 || (nargin == 1 && ischar(AtlasName) && length(AtlasName)>=4 && strcmpi(AtlasName(1:4),'avai'))
    FirstLabel = AtlasList;
    return;
elseif nargin == 1
    if ischar(AtlasName) && any(strcmp(AtlasName,shiSpmAnatLabel))
        PathHere = fileparts(which(['shiSpmTemplate_',AtlasName,'.nii']));
        FirstLabel = fullfile(PathHere,['shiSpmTemplate_',AtlasName,'.nii']);
        AllLabel = fullfile(PathHere,['shiSpmTemplate_',AtlasName,'_Label.txt']);
        return;
    elseif iscell(AtlasName)
        PathHere = cell(size(AtlasName));
        FirstLabel = cell(size(AtlasName));
        AllLabel = cell(size(AtlasName));
        for i = 1:numel(AtlasName)
            PathHere{i} = fileparts(which(['shiSpmTemplate_',AtlasName{i},'.nii']));
            FirstLabel{i} = fullfile(PathHere{i},['shiSpmTemplate_',AtlasName{i},'.nii']);
            AllLabel{i} = fullfile(PathHere{i},['shiSpmTemplate_',AtlasName{i},'_Label.txt']);
        end
        return;
    end
end

if isempty(Roi)
    FirstLabel = cell(0);
    AllLabel = cell(0);
    return;
end

if ~exist('Radius','var') || isempty(Radius)
    Radius = 5;
end

if ~exist('SupressPrint','var') || nargin<3
    SupressPrint = false;
end

if isempty(AtlasName)
    AtlasNum = spm_input('Select an atlas','+1','m',[AtlasList{1},sprintf('|%s',AtlasList{2:end})]);
    AtlasName = AtlasList{AtlasNum};
end

    
% PathHere = fileparts(which(['shiSpmTemplate_',AtlasName,'.img']));
% Template = fullfile(PathHere,['shiSpmTemplate_',AtlasName,'.img']);
% if isempty(PathHere)
    PathHere = fileparts(which(['shiSpmTemplate_',AtlasName,'.nii']));
    Template = fullfile(PathHere,['shiSpmTemplate_',AtlasName,'.nii']);
% end

fid = fopen(fullfile(PathHere,['shiSpmTemplate_',AtlasName,'_Label.txt']));
LabelList=textscan(fid,'%d%s%*[^\n]');
fclose(fid);
LabelNum = LabelList{1};
LabelList = LabelList{2};

% 
% 
% if size(Roi,1)==3 && size(Roi,2)~=3
%     Roi = Roi';
%     warning('Coordinates transposed'); %#ok<WNTAG>
% end
% 
% if ~iscell(Roi)
%     Roi = mat2cell(Roi,ones(size(Roi,1),1),size(Roi,2));
% end
% 
% if ~isnumeric(Roi{1}) && ~ischar(Roi{1})
%     error('unknown ROI format');
% end


Roi = shiSpmRoiFormat(Roi,Radius);

AllLabel = cell(length(Roi),1);

XY=shiSpmRoiXtrVoxelwise(Template,Roi,Radius);

for i = 1:length(AllLabel)

%     Roi1 = Roi{i};
% 
%     xY=shiSpmRoiXtrVoxelwise(Template,Roi1,Radius);
    xY=XY{i};
    sum_n_Voxel = numel(xY);

    xY = xY(~isnan(xY));
    xY = xY(xY~=0);

    value_Voxel = unique(xY);
    n_Voxel = zeros(size(value_Voxel));
    for k = 1:length(n_Voxel)
        n_Voxel(k) = sum(xY==value_Voxel(k));
    end

    [n_Voxel,IND] = sort(n_Voxel,'descend');
    value_Voxel = value_Voxel(IND);
    perc_Voxel = n_Voxel./sum_n_Voxel.*100;

    label_Voxel = cell(size(value_Voxel));
    label_Voxel(:) = LabelList(shiVLookUp(value_Voxel,[LabelNum(:)';1:length(LabelNum)]',2));
    
    if sum_n_Voxel ~= numel(xY)
        n_Voxel(end+1) = sum_n_Voxel - numel(xY); %#ok<*AGROW>
        perc_Voxel(end+1) = 100-sum(perc_Voxel);
        label_Voxel(end+1) = {'unknown'};
    end

    AllLabel{i} = [
        num2cell(perc_Voxel);
        label_Voxel;
        num2cell(n_Voxel);
        ];
    
    if isempty(AllLabel{i})
        AllLabel{i} = {0;'error';0};
    end

    if ~SupressPrint
        disp(Roi{i});
        disp([{'%voxel',AtlasName,'#voxel'};AllLabel{i}']);
    end

end



for i = 1:length(AllLabel)
    FirstLabel{i,1} = AllLabel{i}{2,1};
end

