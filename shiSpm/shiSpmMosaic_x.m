function [MosaicNii,MosaicPng] = shiSpmMosaic_x(ImgName,Gap)

% returns mosaic sagittal nifti images and pictures (.png)
%
% zhenhao shi

if ~exist('Gap','var') || isempty(Gap)
    Gap = 1;
end

cwd=pwd;
ImgName = cellstr(char(ImgName));
MosaicNii = cell(size(ImgName));
MosaicPng = cell(size(ImgName));

for img = 1:length(ImgName)

    xxImg = char(shiFullFileName(ImgName{img}));
    if size(xxImg,1) > 1
        error('do not use wildcard');
    end
    [xxPath,xxName,xxExt] = fileparts(xxImg);
    cd(xxPath);


    V = spm_vol(xxImg);
    Y = spm_read_vols(V);
    Y = Y(1:Gap:end,:,:);

    [x,y,z] = size(Y);

    Xm = ceil(sqrt(x));
    Ym = zeros(Xm*y,Xm*z);  

    for j=1:x
        col = Xm-(floor((j-1)/Xm)+1)+1;
        row = j-floor((j-1)/Xm)*(Xm);
        Ym((row-1)*y+1:row*y,(col-1)*z+1:col*z) = Y(j,:,:);
    end

    clear Vm;
    Vm.fname = ['Mosaic_x_',xxName,xxExt];
    Vm.dim = [size(Ym),1];
    Vm.mat = V.mat;
    Vm.pinfo = V.pinfo;
    Vm.dt = V.dt;
    Vm.n = V.n;
    Vm.descrip = sprintf('Mosaic_x - %s',xxImg);
    spm_write_vol(Vm,Ym);

    MosaicNii{img,1} = fullfile(xxPath,Vm.fname);
    MosaicPng{img,1} = fullfile(xxPath,['Mosaic_x_',xxName,'.png']);

    grey = mat2gray(rot90(Ym));
    grey = kron((grey).*1.7,[1 1;1 1]);
    imwrite(grey,MosaicPng{img,1});

end

cd(cwd);


