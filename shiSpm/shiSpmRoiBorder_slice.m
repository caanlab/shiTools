function ImgRoiBorder = shiSpmRoiBorder_slice(ImgRoi,Coordinate,SliceCut,isInnerBorder,BorderConnectedness,suffix)

% draws ROI border on a particular slice
%
% ImgRoi - ROI image filename
% Coordinate - coordinates at which the slice is cut (e.g. [8,40,0])
% SliceCut - direction of slice cut, 'x' or 'y' or 'z' (e.g. if 'x', then y and z coordinates provided in Coordinate will be ignored)
% isInnerBorder - whether border is just inside ROI (default, true), or just outside ROI (false)
% BorderConnectedness - how to connect border corners (4 or 8) (default = 4)
% suffix - suffix to add to filename (default = '_border')
%

if ~exist('isInnerBorder','var') || isempty(isInnerBorder)
    isInnerBorder = true;
end
if ~exist('BorderConnectedness','var') || isempty(BorderConnectedness)
    BorderConnectedness = 4;
end
if ~exist('suffix','var') || isempty(suffix)
    suffix = '_border';
end

ImgRoi = cellstr(char(ImgRoi));
if length(ImgRoi)>1
    ImgRoiBorder = cell(size(ImgRoi));
    for i = 1:length(ImgRoi)
        ImgRoiBorder{i} = shiSpmRoiBorder_slice(ImgRoi{i},Coordinate,SliceCut,isInnerBorder,BorderConnectedness,suffix);
    end
    return;
end

ImgRoi = char(ImgRoi);

switch BorderConnectedness
    case 4
        Shift = [
            +1 00
            -1 00
            00 +1
            00 -1
            ];
    case 8
        Shift = [
            +1 00
            -1 00
            00 +1
            00 -1
            +1 +1
            -1 +1
            +1 -1
            -1 -1
            ];
end

V = spm_vol(ImgRoi);
[pth,nme,ext] = fileparts(V.fname);
Coordinate = [Coordinate,1];
Index = round(Coordinate/(V.mat)');

Y = spm_read_vols(V);
Y_out = zeros(size(Y));

switch lower(SliceCut)
    case 'x'
        Y_border = shiSpmRoiBorder_slice_shift(squeeze(Y(Index(1),:,:)),Shift,isInnerBorder);
        Y_out(Index(1),:,:) = Y_border;
    case 'y'
        Y_border = shiSpmRoiBorder_slice_shift(squeeze(Y(:,Index(2),:)),Shift,isInnerBorder);
        Y_out(:,Index(2),:) = Y_border;
    case 'z'
        Y_border = shiSpmRoiBorder_slice_shift(squeeze(Y(:,:,Index(3))),Shift,isInnerBorder);
        Y_out(:,:,Index(3)) = Y_border;
end

V.fname = fullfile(pth,[nme,suffix,ext]);
V.descrip = [V.descrip, ' (border)'];

spm_write_vol(V,Y_out);
ImgRoiBorder = V.fname;


function Y_border = shiSpmRoiBorder_slice_shift(Y_plane,Shift,isInnerBorder)

Y_plane_shift = zeros([size(Y_plane),size(Shift,1)]);

for i = 1:size(Shift,1)
    Y_plane_shift(:,:,i) = shiMatShift(Y_plane,0,Shift(i,:));
end

if isInnerBorder
    Y_border = any(Y_plane_shift<0.5,3) & (Y_plane>=0.5);
else
    Y_border = any(Y_plane_shift>=0.5,3) & (Y_plane<0.5);
end
