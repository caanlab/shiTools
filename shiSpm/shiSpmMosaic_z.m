function [MosaicNii,MosaicPng] = shiSpmMosaic_z(ImgName,Gap)

% returns mosaic axial nifti images and pictures (.png)
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
    Y = Y(:,:,1:Gap:end);

    [x,y,z] = size(Y);

    Zm = ceil(sqrt(z));
    Ym = zeros(Zm*x,Zm*y);  

    for j=1:z
        col = Zm-(floor((j-1)/Zm)+1)+1;
        row = j-floor((j-1)/Zm)*(Zm);
        Ym((row-1)*x+1:row*x,(col-1)*y+1:col*y) = Y(:,:,j);
    end

    clear Vm;
    Vm.fname = ['Mosaic_z_',xxName,xxExt];
    Vm.dim = [size(Ym),1];
    Vm.mat = V.mat;
    Vm.pinfo = V.pinfo;
    Vm.dt = V.dt;
    Vm.n = V.n;
    Vm.descrip = sprintf('Mosaic_z - %s',xxImg);
    spm_write_vol(Vm,Ym);

    MosaicNii{img,1} = fullfile(xxPath,Vm.fname);
    MosaicPng{img,1} = fullfile(xxPath,['Mosaic_z_',xxName,'.png']);

    grey = mat2gray(rot90(Ym));
    grey = kron((grey).*1.7,[1 1;1 1]);
    imwrite(grey,MosaicPng{img,1});

end

cd(cwd);


