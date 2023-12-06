function ImgResliced = shiSpmReslice(ImgRaw,InterpolateOrder,Prefix)

% performs reslicing
%
% ImgResliced = shiSpmReslice(ImgRaw)
% ImgResliced = shiSpmReslice(ImgRaw,Prefix)
% ImgResliced = shiSpmReslice(ImgRaw,InterpolateOrder)
% ImgResliced = shiSpmReslice(ImgRaw,InterpolateOrder,Prefix)
%   returns resliced images
%   all resliced images are to be resliced to the first raw image
% 
%   ImgRaw          - raw image file(s)
%   Prefix          - prefix for output resliced image (default = 'Resliced_')
% 
% see spm_mask.m
%
% Zhenhao Shi
% 2022-01-31

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


if ischar(ImgRaw) && size(ImgRaw,1) == 1
    [pt,nme,xt] = fileparts(ImgRaw);
    ImgResliced = fullfile(pt,[Prefix,nme,xt]);
else
    [pt,nme,xt] = shiFileParts(cellstr(char(ImgRaw)));
    ImgResliced = fullfile(pt,shiStrConcat(Prefix,nme,xt));
end

V_Img = spm_vol(P_Img);

m_Img = numel(V_Img);

%-Create headers
%--------------------------------------------------------------------------
V_Out = V_Img;
for i=1:m_Img
    [pth,nm,ext,num] = spm_fileparts(deblank(V_Out(i).fname));
    V_Out(i).fname      = fullfile(pth,[Prefix, nm, ext, num]);
    V_Out(i).descrip    = 'Resliced';
    V_Out(i).mat        = V_Out(1).mat;
    V_Out(i).dim(1:3)   = V_Out(1).dim(1:3);
    if isfield(V_Img(i),'descrip')
        V_Out(i).descrip    = [V_Img(i).descrip,'; Resliced'];
    end
end
V_Out  = spm_create_vol(V_Out);
M   = V_Out(1).mat;
dim = V_Out(1).dim(1:3);

%-Compute resliced images
%--------------------------------------------------------------------------

for j=1:dim(3)

    Mi  = spm_matrix([0 0 j]);

    % Write the images
    for i=1:m_Img
        M1        = M\V_Img(i).mat\Mi;
        img       = spm_slice_vol(V_Img(i),M1,dim(1:2),[InterpolateOrder 0]);
        V_Out(i)  = spm_write_plane(V_Out(i),img,j);
    end

end

