function shiSpmRoiMakeSphere(Coordinate,OutputImg)

% write one spheric ROI in a nifti image with 1x1x1 mm resolution
%
% Coordinate        - n-by-3 or n-by-4 matrix, see shiSpmRoiXtr
% OutputImg         - output image filename(s)
%
% zhenhao shi

if iscell(Coordinate)
    if numel(Coordinate) ~= numel(OutputImg)
        error('number of ROIs must be equal to number of output images');
    end
    for i = 1:numel(Coordinate)
        shiSpmRoiMakeSphere(Coordinate{i},OutputImg{i});
    end
    return;
elseif size(Coordinate,1) > 1
    if size(Coordinate,1) ~= numel(OutputImg)
        error('number of ROIs must be equal to number of output images');
    end
    for i = 1:size(Coordinate,1)
        shiSpmRoiMakeSphere(Coordinate(i,:),OutputImg{i});
    end
    return;
end

V.fname = char(OutputImg);
V.dim = [182 218 182];
V.dt = [2 0];
V.pinfo = [1 1 0]';
V.mat = [
    -1     0     0    91
    0     1     0  -127
    0     0     1   -73
    0     0     0     1
    ];
V.descrip = ['Sphere - ',shiTime];

Y = zeros(V.dim);

[R,C,P]  = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
RCP      = [R(:)';C(:)';P(:)'];
RCP(4,:) = 1;
XYZ      = V.mat(1:3,:) * RCP;

if length(Coordinate) == 3
    Radius = 5;
else
    Radius = Coordinate(4);
end

Distance = sqrt(sum((XYZ'-repmat(Coord,size(XYZ,2),1)).^2,2));
Mask = Distance <= Radius;
Y(Mask) = 1;

spm_write_vol(V,Y);

