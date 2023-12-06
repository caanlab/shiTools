function Img_trim = shiSpmImgTrim(Img)
%
% trims the image to reduce the bounding box so that only non-zero/non-NaN voxels are retained


if iscell(Img) || ~isvector(Img)
    if ischar(Img)
        Img = cellstr(Img);
    end
    Img_trim = cell(size(Img));
    for i = 1:numel(Img)
        Img_trim{i} = shiSpmImgTrim(Img{i});
    end
    return;
end

[pth,nme,ext] = fileparts(Img);
Img_trim = fullfile(pth,[nme,'_trim',ext]);

V = spm_vol(Img);
[Y,XYZ] = spm_read_vols(V);

existNan = any(isnan(Y(:)));

Y0 = false(size(Y));

indx = true(V.dim(1),1);
indy = true(V.dim(2),1);
indz = true(V.dim(3),1);

if ~existNan
    for i = 1:V.dim(1), indx(i) = ~all(Y(i,:,:)==0,'all'); end
    for i = 1:V.dim(2), indy(i) = ~all(Y(:,i,:)==0,'all'); end
    for i = 1:V.dim(3), indz(i) = ~all(Y(:,:,i)==0,'all'); end
else
    for i = 1:V.dim(1), indx(i) = ~all(isnan(Y(i,:,:)),'all'); end
    for i = 1:V.dim(2), indy(i) = ~all(isnan(Y(:,i,:)),'all'); end
    for i = 1:V.dim(3), indz(i) = ~all(isnan(Y(:,:,i)),'all'); end
end

Y0(indx,indy,indz) = true;
Ind_trim = find(Y0);
XYZ_trim = XYZ(:,Ind_trim); %#ok<FNDSB>

Y_trim = Y(indx,indy,indz);

V_trim = V;
V_trim.fname = Img_trim;
V_trim.dim = size(Y_trim);
V_trim.mat(1:3,4) = XYZ_trim(:,1) - sum(V.mat(1:3,1:3),2);

spm_write_vol(V_trim,Y_trim);
