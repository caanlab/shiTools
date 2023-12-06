function shiRestVoxelwiseConn(PathData,PathOutput,ImgMask,doNormalize)

% conducts voxel-wise connectivity analysis and return maps of degrees of various definitions
%
% (masking procedure will be revised)

shiMkdir(PathOutput);

Img = {};
ImgMask = char(ImgMask);

ImgExt = {'.img','.nii','.IMG','.NII'};
for i = 1:length(ImgExt)
    if ~isempty(Img)
        break;
    end
    Img = shiFullFileName(fullfile(PathData,['*',ImgExt{i}]));
end

if isempty(Img)
    error('no data image found');
end
if ~exist(ImgMask,'file')
    error('no mask found');
end

fprintf('\n# of time points = %d\n',length(Img));


%%



% read in 4-D time series

V = spm_vol(char(Img));
[Y0,XYZ] = spm_read_vols(V);


% read in mask

Vm0 = spm_vol(ImgMask);
[Ym0,XYZm0] = spm_read_vols(Vm0);
Ym0(isnan(Ym0)) = 0;
Ym0 = (Ym0 > 0) .* 1;


% reslice mask without writing new file

Xm0_grid = nan(size(Ym0));
Ym0_grid = nan(size(Ym0));
Zm0_grid = nan(size(Ym0));
Xm0_grid(:) = XYZm0(1,:)';
Ym0_grid(:) = XYZm0(2,:)';
Zm0_grid(:) = XYZm0(3,:)';

X_grid = nan(size(Y0(:,:,:,1)));
Y_grid = nan(size(Y0(:,:,:,1)));
Z_grid = nan(size(Y0(:,:,:,1)));
X_grid(:) = XYZ(1,:)';
Y_grid(:) = XYZ(2,:)';
Z_grid(:) = XYZ(3,:)';

Ym = interp3(Ym0_grid,Xm0_grid,Zm0_grid,Ym0,Y_grid,X_grid,Z_grid,'nearest',0);
n_Vox = sum(Ym(:));
fprintf('# of in-mask voxels = %d   ',n_Vox);


% apply mask

Y = Y0.*repmat(Ym,[1,1,1,size(Y0,4)]);


% pick up voxels in mask

vecYm = Ym(:)';

vecY = nan(size(Y,4),size(Y,1)*size(Y,2)*size(Y,3));
for t = 1:size(Y,4)
    tmp = Y(:,:,:,t);
    vecY(t,:) = tmp(:);
end

indInMask = find(vecYm>0);


% calculate metrics based on corr coef

Y_DegreePos = nan(size(Ym));
Y_DegreeNeg = nan(size(Ym));
Y_DegreeAbs = nan(size(Ym));
Y_PowerPos = nan(size(Ym));
Y_PowerNeg = nan(size(Ym));
Y_PowerAll = nan(size(Ym));

fName_DegreePos = 'vxw_0001_DegreePos.img';
fName_DegreeNeg = 'vxw_0002_DegreeNeg.img';
fName_DegreeAbs = 'vxw_0003_DegreeAbs.img';
fName_PowerPos  = 'vxw_0004_PowerPos.img';
fName_PowerNeg  = 'vxw_0005_PowerNeg.img';
fName_PowerAll  = 'vxw_0006_PowerAll.img';

Vtemplate.dt = [64 0];
Vtemplate.dim = size(Ym);
Vtemplate.pinfo = [1 0 0]';
Vtemplate.mat = V(1).mat;
Vtemplate.descrip = 'shiRestVoxelwiseConn';

V_DegreePos = Vtemplate;  
V_DegreeNeg = Vtemplate;
V_DegreeAbs = Vtemplate;
V_PowerPos  = Vtemplate;
V_PowerNeg  = Vtemplate;
V_PowerAll  = Vtemplate;

V_DegreePos.fname = fullfile(PathOutput,fName_DegreePos);
V_DegreeNeg.fname = fullfile(PathOutput,fName_DegreeNeg);
V_DegreeAbs.fname = fullfile(PathOutput,fName_DegreeAbs);
V_PowerPos.fname  = fullfile(PathOutput,fName_PowerPos );
V_PowerNeg.fname  = fullfile(PathOutput,fName_PowerNeg );
V_PowerAll.fname  = fullfile(PathOutput,fName_PowerAll );

fprintf('00.00%%');
cnt = 0;

for i = 1:numel(Ym)

    if Ym(i) == 0
        continue;
    end
    
    cnt = cnt + 1;
    fprintf('\b\b\b\b\b\b%5.2f%%',cnt/n_Vox*100);

    xR = corr(vecY(:,setdiff(indInMask,i)),vecY(:,i));
    
    xZ = atanh(xR);
    xZ(isinf(xZ)) = NaN;
    xR2 = xR.^2;

    xZ_Pos = xZ(xR>0);
    xZ_Neg = xZ(xR<0);
    xZ_Abs = abs(xZ);
    xR2_Pos = xR2(xR>0);
    xR2_Neg = xR2(xR<0);
    
    Y_DegreePos(i) = nanmean(xZ_Pos);
    Y_DegreeNeg(i) = nanmean(xZ_Neg);
    Y_DegreeAbs(i) = nanmean(xZ_Abs);
    Y_PowerPos(i) = nanmean(xR2_Pos);
    Y_PowerNeg(i) = nanmean(xR2_Neg);
    Y_PowerAll(i) = nanmean(xR2);

end


% finish

fprintf('\n# of voxels in output images\n');
    fprintf('   (DegreePos) = %d\n',nansum(Y_DegreePos(:)<Inf));
    fprintf('   (DegreeNeg) = %d\n',nansum(Y_DegreeNeg(:)<Inf));
    fprintf('   (DegreeAbs) = %d\n',nansum(Y_DegreeAbs(:)<Inf));
    fprintf('   (PowerPos)  = %d\n',nansum(Y_PowerPos(:)<Inf));
    fprintf('   (PowerNeg)  = %d\n',nansum(Y_PowerNeg(:)<Inf));
    fprintf('   (PowerAll)  = %d\n',nansum(Y_PowerAll(:)<Inf));

spm_write_vol(V_DegreePos,Y_DegreePos);
spm_write_vol(V_DegreeNeg,Y_DegreeNeg);
spm_write_vol(V_DegreeAbs,Y_DegreeAbs);
spm_write_vol(V_PowerPos ,Y_PowerPos );
spm_write_vol(V_PowerNeg ,Y_PowerNeg );
spm_write_vol(V_PowerAll ,Y_PowerAll );


% normalize

if ~doNormalize
    return;
end

Vxw_DegreePos_div = fullfile(PathOutput,['vxw_0001_DegreePos_div','.img']);
Vxw_DegreeNeg_div = fullfile(PathOutput,['vxw_0002_DegreeNeg_div','.img']);
Vxw_DegreeAbs_div = fullfile(PathOutput,['vxw_0003_DegreeAbs_div','.img']);
Vxw_PowerPos_div  = fullfile(PathOutput,['vxw_0004_PowerPos_div' ,'.img']);
Vxw_PowerNeg_div  = fullfile(PathOutput,['vxw_0005_PowerNeg_div' ,'.img']);
Vxw_PowerAll_div  = fullfile(PathOutput,['vxw_0006_PowerAll_div' ,'.img']);

Vxw_DegreePos_z = fullfile(PathOutput,['vxw_0001_DegreePos_z','.img']);
Vxw_DegreeNeg_z = fullfile(PathOutput,['vxw_0002_DegreeNeg_z','.img']);
Vxw_DegreeAbs_z = fullfile(PathOutput,['vxw_0003_DegreeAbs_z','.img']);
Vxw_PowerPos_z  = fullfile(PathOutput,['vxw_0004_PowerPos_z' ,'.img']);
Vxw_PowerNeg_z  = fullfile(PathOutput,['vxw_0005_PowerNeg_z' ,'.img']);
Vxw_PowerAll_z  = fullfile(PathOutput,['vxw_0006_PowerAll_z' ,'.img']);

MEAN_DegreePos = num2str(nanmean(Y_DegreePos(:)));
MEAN_DegreeNeg = num2str(nanmean(Y_DegreeNeg(:)));
MEAN_DegreeAbs = num2str(nanmean(Y_DegreeAbs(:)));
MEAN_PowerPos  = num2str(nanmean(Y_PowerPos(:)));
MEAN_PowerNeg  = num2str(nanmean(Y_PowerNeg(:)));
MEAN_PowerAll  = num2str(nanmean(Y_PowerAll(:)));

SD_DegreePos = num2str(nanstd(Y_DegreePos(:)));
SD_DegreeNeg = num2str(nanstd(Y_DegreeNeg(:)));
SD_DegreeAbs = num2str(nanstd(Y_DegreeAbs(:)));
SD_PowerPos  = num2str(nanstd(Y_PowerPos(:)));
SD_PowerNeg  = num2str(nanstd(Y_PowerNeg(:)));
SD_PowerAll  = num2str(nanstd(Y_PowerAll(:)));


shiSpmImgCalc0(fullfile(PathOutput,fName_DegreePos),['i1./',MEAN_DegreePos],Vxw_DegreePos_div);
shiSpmImgCalc0(fullfile(PathOutput,fName_DegreeNeg),['i1./',MEAN_DegreeNeg],Vxw_DegreeNeg_div);
shiSpmImgCalc0(fullfile(PathOutput,fName_DegreeAbs),['i1./',MEAN_DegreeAbs],Vxw_DegreeAbs_div);
shiSpmImgCalc0(fullfile(PathOutput,fName_PowerPos ),['i1./',MEAN_PowerPos ],Vxw_PowerPos_div );
shiSpmImgCalc0(fullfile(PathOutput,fName_PowerNeg ),['i1./',MEAN_PowerNeg ],Vxw_PowerNeg_div );
shiSpmImgCalc0(fullfile(PathOutput,fName_PowerAll ),['i1./',MEAN_PowerAll ],Vxw_PowerAll_div );

shiSpmImgCalc0(fullfile(PathOutput,fName_DegreePos),['(i1-',MEAN_DegreePos,')./',SD_DegreePos],Vxw_DegreePos_z);
shiSpmImgCalc0(fullfile(PathOutput,fName_DegreeNeg),['(i1-',MEAN_DegreeNeg,')./',SD_DegreeNeg],Vxw_DegreeNeg_z);
shiSpmImgCalc0(fullfile(PathOutput,fName_DegreeAbs),['(i1-',MEAN_DegreeAbs,')./',SD_DegreeAbs],Vxw_DegreeAbs_z);
shiSpmImgCalc0(fullfile(PathOutput,fName_PowerPos ),['(i1-',MEAN_PowerPos,')./',SD_PowerPos ],Vxw_PowerPos_z );
shiSpmImgCalc0(fullfile(PathOutput,fName_PowerNeg ),['(i1-',MEAN_PowerNeg,')./',SD_PowerNeg ],Vxw_PowerNeg_z );
shiSpmImgCalc0(fullfile(PathOutput,fName_PowerAll ),['(i1-',MEAN_PowerAll,')./',SD_PowerAll ],Vxw_PowerAll_z );

shiSpmMask(Vxw_DegreePos_z,ImgMask,'inclusive','');
shiSpmMask(Vxw_DegreeNeg_z,ImgMask,'inclusive','');
shiSpmMask(Vxw_DegreeAbs_z,ImgMask,'inclusive','');
shiSpmMask(Vxw_PowerPos_z ,ImgMask,'inclusive','');
shiSpmMask(Vxw_PowerNeg_z ,ImgMask,'inclusive','');
shiSpmMask(Vxw_PowerAll_z ,ImgMask,'inclusive','');
    
