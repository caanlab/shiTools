function OutputImgName = shiSpmRoiDraw(Roi,OutputImgName,Space)

% writes Roi into binary (0 vs 1) nifti file(s)
%
% Roi can be a binary .img/.nii image, or MNI coordinates ([x,y,z]) or MNI coordinates+radius ([x,y,z,radius])
% 
% OutputImgName = shiSpmRoiDraw(Roi,OutputImgName)
%
%    ###########
% by Zhenhao Shi @ 2018-8-7
% 

thresh = 0;

[Roi,RoiFormat] = shiSpmRoiFormat(Roi);

if ~exist('OutputImgName','var') || isempty(OutputImgName)
    OutputImgName = shiStrConcat('Roi_',1:length(Roi),'.nii');
end
OutputImgName = cellstr(char(OutputImgName));

if length(Roi) ~= length(OutputImgName)
    error('number of Rois must be equal to number of output image file names');
end

if ~exist('Space','var')
    Space = '';
end


for r = 1:length(Roi)

    if isequal(RoiFormat{r},'Mask')

        [Mat_Img,Dim_Img] = xxSpace2V(Space);
        msk = false(Dim_Img);

        %% read ROI file

        File_Roi = Roi{r};
        V_Roi = spm_vol(File_Roi);
        Mat_Roi = V_Roi.mat;

        %% compute ROI mask slice-by-slice in Image space

        for j = 1:Dim_Img(3)

            x_Slice_Img = spm_matrix([0 0 j]); % matrix, for a slice of Image

            Mat_Roi_to_ImgSlice = Mat_Img\Mat_Roi\x_Slice_Img; % matrix, from roi to a slice of Image
            img_slice = spm_slice_vol(V_Roi,Mat_Roi_to_ImgSlice,Dim_Img(1:2),[0 NaN]); % a slice of ROI in Image space

            msk(:,:,j) = img_slice > thresh;

        end

    elseif isequal(RoiFormat{r},'Atlas')

        if isempty(Space)
            xSpace = Roi{r}.atlas;
        else
            xSpace = Space;
        end
        [Mat_Img,Dim_Img] = xxSpace2V(xSpace);
        msk = false(Dim_Img);

        %% read ROI file

        File_Roi = Roi{r}.atlas;
        V_Roi = spm_vol(File_Roi);
        Mat_Roi = V_Roi.mat;

        %% compute ROI mask slice-by-slice in Image space

        for j = 1:Dim_Img(3)

            x_Slice_Img = spm_matrix([0 0 j]); % matrix, for a slice of Image

            Mat_Roi_to_ImgSlice = Mat_Img\Mat_Roi\x_Slice_Img; % matrix, from roi to a slice of Image
            img_slice = spm_slice_vol(V_Roi,Mat_Roi_to_ImgSlice,Dim_Img(1:2),[0 NaN]); % a slice of ROI in Image space

            msk(:,:,j) = img_slice == Roi{r}.value;

        end

    elseif isequal(RoiFormat{r},'Sphere')

        [Mat_Img,Dim_Img,XYZ_Img] = xxSpace2V(Space);
        msk = false(Dim_Img);

        %% make sphere

        distance = sqrt(sum((XYZ_Img'-repmat(Roi{r}(1:3),size(XYZ_Img,2),1)).^2,2));
        msk_index = distance<=Roi{r}(4);
        if ~any(msk_index)
            msk_index = distance == min(distance);
        end

        msk(msk_index) = true;

    else

        error('Roi # %d: unknown Roi format',r);

    end

    xV = struct('fname',OutputImgName{r},'dim',Dim_Img,'dt',[2,0],'mat',Mat_Img);
    spm_write_vol(xV,msk);

end


function [Mat_Img,Dim_Img,XYZ_Img] = xxSpace2V(Space)

if isempty(Space)
    Space = which('shi_defaultUnderlay.nii');
end
if exist(Space,'file')
    V_Img = spm_vol(Space);
elseif isstruct(Space)
    V_Img = Space;
else
    error('Space must be a nifti image or "V" variable from spm_vol');
end

Mat_Img = V_Img(1).mat;
Dim_Img = V_Img(1).dim(1:3);

if nargout==3
    [R,C,P] = ndgrid(1:Dim_Img(1),1:Dim_Img(2),1:Dim_Img(3)); % see spm_read_vols
    RCP     = [R(:)';C(:)';P(:)';ones(1,numel(R))];
    XYZ_Img   = Mat_Img(1:3,:)*RCP;
end