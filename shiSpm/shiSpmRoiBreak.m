function ImgOut = shiSpmRoiBreak(Img,Conn)

% finds disconnected components in a binary/mask image and save them as individual images

if ~exist('Conn','var') || isempty(Conn) || (Conn-6)*(Conn-18)*(Conn-26) ~= 0
    Conn = 18;
end

if iscellstr(Img)
    ImgOut = cell(size(Img));
    for i = 1:numel(Img)
        ImgOut{i} = shiSpmRoiBreak(Img{i},Conn);
    end
elseif ~ischar(Img)
    error('Img must be char or cellstr');
end

[pt,nm,xt] = fileparts(Img);
pt2 = pt;

if isempty(pt)
    pt2 = pwd;
end

V = spm_vol(Img);
Y = spm_read_vols(V);

[YL,nCluster] = spm_bwlabel(Y,Conn);

ImgOut = shiStrConcat(pt2,filesep,nm,'_comp',1:nCluster,xt);

for i = 1:nCluster
    xV = V;
    xV.fname = ImgOut{i};
    xV.descrip = sprintf('%s - (component %d/%d)',xV.descrip,i,nCluster);
    xY = YL == i;
    spm_write_vol(xV,xY);
end


