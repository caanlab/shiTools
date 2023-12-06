function [V_ImgOut,Y_ImgOut,XYZ] = shiSpmResliceRead(ImgRaw,InterpolateOrder,Prefix)

% performs reslicing without saving resliced image files
%
% [V_ImgOut,Y_ImgOut,XYZ] = shiSpmResliceRead(ImgRaw)
% [V_ImgOut,Y_ImgOut,XYZ] = shiSpmResliceRead(ImgRaw,Prefix)
% [V_ImgOut,Y_ImgOut,XYZ] = shiSpmResliceRead(ImgRaw,InterpolateOrder)
% [V_ImgOut,Y_ImgOut,XYZ] = shiSpmResliceRead(ImgRaw,InterpolateOrder,Prefix)
%   returns resliced V, Y, and XYZ variables (see spm_vol.m)
%   all images are to be resliced to the first raw image
% 
%   ImgRaw          - raw image file(s)
%   Prefix          - prefix for output resliced image (default = 'Resliced_')
%   V_ImgOut/Y_ImgOut           - V/Y variable for resliced image output
%   XYZ                         - coordinates
%
% see spm_mask.m
%
% Zhenhao Shi
% 2018-04-11


if ~nargin
    [P_Img, sts] = spm_select(Inf,'image','Images to reslice');
    if ~sts, return; end
else
    P_Img = char(cellstr(ImgRaw));
end

if exist('InterpolateOrder','var') && ischar(InterpolateOrder) && ~exist('Prefix','var')
    Prefix = InterpolateOrder;
    InterpolateOrder = 1;
elseif ~exist('InterpolateOrder','var') || isempty(InterpolateOrder)
    InterpolateOrder = 1;
end

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'Resliced_';
end



V_Img = spm_vol(P_Img);

m_Img = numel(V_Img);

%-Create headers
%--------------------------------------------------------------------------
V_ImgOut = V_Img;
for i=1:m_Img
    [pth,nm,ext,num] = spm_fileparts(deblank(V_ImgOut(i).fname));
    V_ImgOut(i).fname      = fullfile(pth,[Prefix, nm, ext, num]);
    V_ImgOut(i).descrip    = 'Resliced';
    V_ImgOut(i).mat        = V_ImgOut(1).mat;
    V_ImgOut(i).dim(1:3)   = V_ImgOut(1).dim(1:3);
    if isfield(V_Img(i),'descrip')
        V_ImgOut(i).descrip    = [V_Img(i).descrip,'; Resliced'];
    end
end

% V_ImgOut  = spm_create_vol(V_ImgOut);
M   = V_ImgOut(1).mat;
dim = V_ImgOut(1).dim(1:3);
[Y_Out1,XYZ] = spm_read_vols(V_Img(1));
Y_ImgOut = nan([size(Y_Out1),m_Img]);

%-Compute resliced images
%--------------------------------------------------------------------------

for j=1:dim(3)

    Mi  = spm_matrix([0 0 j]);
    
    % Write the images
    for i=1:m_Img
        M1        = M\V_Img(i).mat\Mi;
        img       = spm_slice_vol(V_Img(i),M1,dim(1:2),[InterpolateOrder 0]);
        Y_ImgOut(:,:,j,i) = img;
    end

end

