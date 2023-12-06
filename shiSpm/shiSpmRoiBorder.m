function ImgRoiBorder = shiSpmRoiBorder(ImgRoi,isInnerBorder,BorderConnectedness,suffix)

% draws 3-dimensional ROI border
%
% ImgRoi - ROI image filename
% isInnerBorder - whether border is just inside ROI (default), or just outside ROI
% BorderConnectedness - how to connect border corners (6, 18, or 26) (default = 6)
% suffix - suffix to add to filename (default = '_border')
%

if ~exist('isInnerBorder','var') || isempty(isInnerBorder)
    isInnerBorder = true;
end
if ~exist('BorderConnectedness','var') || isempty(BorderConnectedness)
    BorderConnectedness = 6;
end
if ~exist('suffix','var') || isempty(suffix)
    suffix = '_border';
end

ImgRoi = cellstr(char(ImgRoi));
if length(ImgRoi)>1
    ImgRoiBorder = cell(size(ImgRoi));
    for i = 1:length(ImgRoi)
        ImgRoiBorder{i} = shiSpmRoiBorder(ImgRoi{i},isInnerBorder,BorderConnectedness,suffix);
    end
    return;
end

ImgRoi = char(ImgRoi);

switch BorderConnectedness
    case 6
        Shift = [
            +1 00 00
            -1 00 00
            00 +1 00
            00 -1 00
            00 00 +1
            00 00 -1
            ];
    case 18
        Shift = [
            +1 00 00
            -1 00 00
            00 +1 00
            00 -1 00
            00 00 +1
            00 00 -1

            +1 +1 00
            -1 -1 00
            +1 -1 00
            -1 +1 00

            00 +1 +1
            00 -1 -1
            00 +1 -1
            00 -1 +1

            +1 00 +1
            -1 00 -1
            -1 00 +1
            +1 00 -1
            ];
    case 26
        Shift = [
            +1 00 00
            -1 00 00
            00 +1 00
            00 -1 00
            00 00 +1
            00 00 -1

            +1 +1 00
            -1 -1 00
            +1 -1 00
            -1 +1 00

            00 +1 +1
            00 -1 -1
            00 +1 -1
            00 -1 +1

            +1 00 +1
            -1 00 -1
            -1 00 +1
            +1 00 -1
            
            +1 +1 +1
            +1 +1 -1
            +1 -1 +1
            +1 -1 -1
            -1 +1 +1
            -1 +1 -1
            -1 -1 +1
            -1 -1 -1
            ];
end

V = spm_vol(ImgRoi);
[pth,nme,ext] = fileparts(V.fname);
Y = spm_read_vols(V);
Y_shift = zeros([size(Y),size(Shift,1)]);

for i = 1:size(Shift,1)
    Y_shift(:,:,:,i) = shiMatShift(Y,0,Shift(i,:));
end

if isInnerBorder
    Y_border = any(Y_shift<0.5,4) & (Y>=0.5);
else
    Y_border = any(Y_shift>=0.5,4) & (Y<0.5);
end

V.fname = fullfile(pth,[nme,suffix,ext]);
V.descrip = [V.descrip, ' (border)'];

spm_write_vol(V,Y_border);
ImgRoiBorder = V.fname;
