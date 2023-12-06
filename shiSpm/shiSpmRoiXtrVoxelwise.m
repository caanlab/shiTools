function [X,Coord,Roi,RoiFormat] = shiSpmRoiXtrVoxelwise(Img,Roi,Radius)

% returns values of each voxel in ROIs
%
% see shiSpmRoiXtr for "Roi" variable format
% shiSpmRoiXtrVoxelwise(Img,Roi)
% shiSpmRoiXtrVoxelwise(Img,Roi,Radius) (Radius does NOT overwrite Roi)
% X: cell array, one cell for each Roi. Each cell is a n_Img-by-n_Voxel matrix
%
% zhenhao shi 2018-7-3

thresh = 0;

if ~exist('Radius','var') || isempty(Radius) || ~isnumeric(Radius) || ~(Radius>=0)
    Radius = 5;
end

%% check Img 

File_Img = char(Img);
for i = size(File_Img,1)
    if ~exist(deblank(File_Img(i,:)),'file')
        error('cannot find %s',File_Img(i,:))
    end
end

if ~exist('Roi','var') || isempty(Roi)
    V_Img = spm_vol(File_Img);
    [Y_Img,XYZ_Img] = spm_read_vols(V_Img);
    X = {reshape(Y_Img,[prod(V_Img(1).dim(1:3)),size(File_Img,1)])'};
    Coord = {XYZ_Img'};
    Roi = '';
    RoiFormat = '';
    return;
end
    


%%
[Roi,RoiFormat] = shiSpmRoiFormat(Roi,Radius);


%%
% %% convert vector and string to cell
% 
% if isnumeric(Roi)
%     Roi = mat2cell(Roi,ones(size(Roi,1),1),size(Roi,2));
% elseif ischar(Roi)
%     Roi = cellstr(Roi);
% end
% 
% 
% %% check each ROI
% 
% RoiFormat = cell(length(Roi),1);
% for r = 1:length(Roi)
%     if isnumeric(Roi{r})
%         if ~isequal(size(Roi{r}),[1,3]) && ~isequal(size(Roi{r}),[1,4])
%             error('Roi coordinate information must be [x,y,z,radius] or [x,y,z] (default radius = 5)');
%         elseif size(Roi{r},2)==3
%             Roi{r} = [Roi{r},Radius];
%         elseif ~(Roi{r}(4)>=0)
%             error('Roi radius must be >= 0 (when =0, extract nearest voxel)');
%         end
%         RoiFormat{r} = 'MNI';
%     elseif ischar(Roi{r})
%         if ~exist(Roi{r},'file')
%             error('cannot find %s',Roi{r});
%         end
%         RoiFormat{r} = 'SingleLabel';
%     else
%         error('unknown ROI format');
%     end
% end


%%
%%

V_Img = spm_vol(File_Img);
[Y_Img,XYZ_Img] = spm_read_vols(V_Img);

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

X = cell(length(Roi),1);
Coord = cell(length(Roi),1);

for r = 1:length(Roi)

    if isequal(RoiFormat{r},'Mask')


        %% read ROI file
        
        File_Roi = Roi{r};
        V_Roi = spm_vol(File_Roi);
        Mat_Roi = V_Roi.mat;


        %% compute ROI mask slice-by-slice in Image space

        msk = true(Dim_Img);

        for j = 1:Dim_Img(3)

            x_Slice_Img = spm_matrix([0 0 j]); % matrix, for a slice of Image

            Mat_Roi_to_ImgSlice = Mat_Img\Mat_Roi\x_Slice_Img; % matrix, from roi to a slice of Image
            img_slice = spm_slice_vol(V_Roi,Mat_Roi_to_ImgSlice,Dim_Img(1:2),[0 NaN]); % a slice of ROI in Image space

            msk(:,:,j) = img_slice > thresh;

        end


        %% apply mask

        msk = reshape(msk,prod(Dim_Img),1)';

        X{r} = Y_Img(:,msk);
        Coord{r} = XYZ_Img(:,msk)';

    elseif isequal(RoiFormat{r},'Atlas')


        %% read ROI file
        
        File_Roi = Roi{r}.atlas;
        V_Roi = spm_vol(File_Roi);
        Mat_Roi = V_Roi.mat;


        %% compute ROI mask slice-by-slice in Image space

        msk = true(Dim_Img);

        for j = 1:Dim_Img(3)

            x_Slice_Img = spm_matrix([0 0 j]); % matrix, for a slice of Image

            Mat_Roi_to_ImgSlice = Mat_Img\Mat_Roi\x_Slice_Img; % matrix, from roi to a slice of Image
            img_slice = spm_slice_vol(V_Roi,Mat_Roi_to_ImgSlice,Dim_Img(1:2),[0 NaN]); % a slice of ROI in Image space

            msk(:,:,j) = img_slice == Roi{r}.value;

        end


        %% apply mask

        msk = reshape(msk,prod(Dim_Img),1)';

        X{r} = Y_Img(:,msk);
        Coord{r} = XYZ_Img(:,msk)';

    elseif isequal(RoiFormat{r},'Sphere')


        %% make sphere
        
        distance = sqrt(sum((XYZ_Img'-repmat(Roi{r}(1:3),size(XYZ_Img,2),1)).^2,2));
        msk = distance<=Roi{r}(4);
        if ~any(msk)
            msk = distance == min(distance);
        end
        X{r} = Y_Img(:,msk);
        Coord{r} = XYZ_Img(:,msk)';


    else
        
        error('Roi # %d: unknown ROI format',r);
        
    end

end
