function [MosaicNii,MosaicPng] = shiSpmMosaic_y(ImgName,Gap)

% returns mosaic coronal nifti images and pictures (.png)
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
    Y = Y(:,1:Gap:end,:);

    [x,y,z] = size(Y);

    Wm = ceil(sqrt(y));
    Ym = zeros(Wm*x,Wm*z);  

    for j=1:y
        col = Wm-(floor((j-1)/Wm)+1)+1;
        row = j-floor((j-1)/Wm)*(Wm);
        Ym((col-1)*x+1:col*x,(row-1)*z+1:row*z) = squeeze(Y(:,j,:));
    end

    clear Vm;
    Vm.fname = ['Mosaic_y_',xxName,xxExt];
    Vm.dim = [size(Ym),1];
    Vm.mat = V.mat;
    Vm.pinfo = V.pinfo;
    Vm.dt = V.dt;
    Vm.n = V.n;
    Vm.descrip = sprintf('Mosaic_y - %s',xxImg);
    spm_write_vol(Vm,Ym);

    MosaicNii{img,1} = fullfile(xxPath,Vm.fname);
    MosaicPng{img,1} = fullfile(xxPath,['Mosaic_y_',xxName,'.png']);

    grey = mat2gray(rot90(Ym));
    grey = kron((grey).*1.7,[1 1;1 1]);
    imwrite(grey,MosaicPng{img,1});

end

cd(cwd);


