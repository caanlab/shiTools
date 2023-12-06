function shiSpmRoiMakeMultisphere(NodeNumber,NodeCoordinate,OutputImg)

% writes multiple spheric ROIs in a single nifti image with 1x1x1 mm resolution
%
% NodeNumber        - one value for each Roi, written into the output image
% NodeCoordinate    - n-by-3 or n-by-4 matrix, see shiSpmRoiXtr
% OutputImg         - output image filename
%
% zhenhao shi

V.fname = char(OutputImg);
V.dim = [182 218 182];
V.dt = [8 0];
V.pinfo = [1 1 0]';
V.mat = [
    -1     0     0    91
    0     1     0  -127
    0     0     1   -73
    0     0     0     1
    ];
V.descrip = ['multisphere - ',shiTime];

Y = zeros(V.dim);

[R,C,P]  = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
RCP      = [R(:)';C(:)';P(:)'];
RCP(4,:) = 1;
XYZ      = V.mat(1:3,:) * RCP;


if size(NodeCoordinate,2) == 3
    NodeCoordinate(:,4) = 5;
end

for i = 1:length(NodeNumber)
    fprintf('.');
    xCoord = NodeCoordinate(i,1:3);
    xRadius = NodeCoordinate(i,4);
    xDistance = sqrt(sum((XYZ'-repmat(xCoord,size(XYZ,2),1)).^2,2));
    xMask = xDistance <= xRadius;
    Y(xMask) = NodeNumber(i);
end

spm_write_vol(V,Y);

