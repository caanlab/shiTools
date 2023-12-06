function shiSpmAlphaSim_RegPerm(Image, X, VarOfInterest, WorkingDir, MaskImg, OutTxtPrefix, ConnectDef, Voxel_P, n_Iteration)

% performs permutation-based Monte-Carlo simulation using the "Eklund" approach
%
% (masking procedure will be revised)


[Image,X] = shi_deNaN(Image,X);
nObs = size(X,1);
nVar = size(X,2);

if any(VarOfInterest>nVar) || any(VarOfInterest<=0)
    error('');
end

DF = nObs-nVar-1;
Voxel_T = tinv(1-Voxel_P,DF);

MaskImg_Resliced = char(shiSpmReslice(MaskImg,Image{1},WorkingDir,0,'perm_RelicedMask_'));
Y_MaskImg_Resliced = shiNiftiRead(MaskImg_Resliced) > 0;


for i = VarOfInterest

    Dir_root = shiMkdir(fullfile(WorkingDir,['perm_Vec',num2str(i,'%2d')]));
    X_reduced = X(:,setdiff(1:nVar,i)); % X

    [~,~,residImage_reduced] = shiSpmAlphaSim_RegPerm_regress(Dir_root,Image,X_reduced);
    fittedImage_reduced = shiSpmAlphaSim_RegPerm_calc_pairwise(Dir_root,Image,residImage_reduced,'i1-i2','fitted_reduced_');

    MaxClusterSize_Pos = nan(n_Iteration,numel(Voxel_T));
    MaxClusterSize_Neg = nan(n_Iteration,numel(Voxel_T));
    MaxT_Pos = nan(n_Iteration,1);
    MaxT_Neg = nan(n_Iteration,1);

    for j = 1:n_Iteration
        
        shiDisp({['Var: ',num2str(i)],['Iter: ',num2str(j)]});

        Dir_perm = shiMkdir(fullfile(Dir_root,['Perm_',num2str(j,'%5d')]));
        residImage_reduced_perm = residImage_reduced(randperm(nObs));

        Image_perm = shiSpmAlphaSim_RegPerm_calc_pairwise(Dir_perm,fittedImage_reduced,residImage_reduced_perm,'i1+i2','perm_');

        [~,tImage_perm,~] = shiSpmAlphaSim_RegPerm_regress(Dir_perm,Image_perm,X);

        Y_TImg = shiNiftiRead(tImage_perm{i});
        Y_TImg(~Y_MaskImg_Resliced) = 0;

        tmp_Pos = nan(1,numel(Voxel_T));
        tmp_Neg = nan(1,numel(Voxel_T));
        for p = 1:numel(Voxel_T)
            Y_TImg_masked_Pos = ( Y_TImg > Voxel_T(p) ) + 0;
            Y_TImg_masked_Neg = ( Y_TImg < -Voxel_T(p) ) + 0;
            tmp_Pos(1,p) = shiSpmAlphaSim_RegPerm_MaxCluster(Y_TImg_masked_Pos,ConnectDef);
            tmp_Neg(1,p) = shiSpmAlphaSim_RegPerm_MaxCluster(Y_TImg_masked_Neg,ConnectDef);
        end
        MaxClusterSize_Pos(j,:) = tmp_Pos;
        MaxClusterSize_Neg(j,:) = tmp_Neg;
        MaxT_Pos(j) = nanmax(Y_TImg(:));
        MaxT_Neg(j) = -nanmin(Y_TImg(:));

        rmdir(Dir_perm,'s');

    end


    cd(WorkingDir);

    %%
    MaxClusterSize_sorted_Pos = MaxClusterSize_Pos;
    MaxClusterSize_sorted_Neg = MaxClusterSize_Neg;
    fid = fopen([OutTxtPrefix,num2str(i,'%2d'),'.txt'],'a+');
    fprintf(fid,'\n\n%s\n%s\n\n',repmat('#',1,19),shiTime(0));
    for p = 1:numel(Voxel_T)
        MaxClusterSize_sorted_Pos(:,p) = sort(MaxClusterSize_sorted_Pos(:,p),'descend');
        fprintf(fid,'p_voxel_Pos<%g\t',Voxel_P(p));
    end
    for p = 1:numel(Voxel_T)
        MaxClusterSize_sorted_Neg(:,p) = sort(MaxClusterSize_sorted_Neg(:,p),'descend');
        fprintf(fid,'p_voxel_Neg<%g\t',Voxel_P(p));
    end
    fprintf(fid,'Tmax_Pos\tTmax_Neg\tp_Cluster\n\n');
    fclose(fid);
    MaxClusterSize_Pos_sorted = [MaxClusterSize_sorted_Pos';MaxClusterSize_sorted_Neg';sort(MaxT_Pos,'descend')';sort(MaxT_Neg,'descend')';(1:n_Iteration)/n_Iteration]';
    dlmwrite([OutTxtPrefix,num2str(i,'%2d'),'.txt'],MaxClusterSize_Pos_sorted,'-append','delimiter','\t');


end




function [betaImage,tImage,residImage] = shiSpmAlphaSim_RegPerm_regress(Dir,Img,Cov)

CovName = cell(size(Cov,2));
for j = 1:size(Cov,2)
    CovName{j} = ['Vector',num2str(j)];
end

matlabbatch{1}.spm.stats.factorial_design.dir = {Dir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = Img;


for j = 1:size(Cov,2)
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).c = Cov(:,j);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).cname = CovName{j};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).iCC = 5;
end
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).c = ones(length(Img),1);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).cname = 'Constant';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).iCC = 5;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');

for j = 1:size(Cov,2)
    matlabbatch{3}.spm.stats.con.consess{j}.tcon.name = [num2str(j,'%.2d'),'_Pos_',CovName{j}];
    matlabbatch{3}.spm.stats.con.consess{j}.tcon.convec = [zeros(1,j-1),1,zeros(1,size(Cov,2)-j),0];
    matlabbatch{3}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
end

matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('serial',matlabbatch);

betaImage = shiFullFileName(fullfile(Dir,'beta_*.nii'));
tImage = shiFullFileName(fullfile(Dir,'spmT_*.nii'));
residImage = shiFullFileName(fullfile(Dir,'Res_*.nii'));

if isempty(betaImage)
    betaImage = shiFullFileName(fullfile(Dir,'beta_*.img'));
    tImage = shiFullFileName(fullfile(Dir,'spmT_*.img'));
    residImage = shiFullFileName(fullfile(Dir,'Res_*.img'));
end



function outImg = shiSpmAlphaSim_RegPerm_calc_pairwise(Dir,Img1,Img2,Func,Prefix)
if numel(Img1) ~= numel(Img2)
    error('');
end
outImg = cell(size(Img1));
for i = 1:numel(Img1)
    outImg{i} = fullfile(Dir,[Prefix,num2str(i,'%2d'),'.nii']);
    shiSpmImgCalc0({Img1{i};Img2{i}},Func,outImg{i});
end




function MaxCluster = shiSpmAlphaSim_RegPerm_MaxCluster(Y_TImg_masked,ConnectDef)

MaxCluster = 0;

[L,num] = spm_bwlabel(Y_TImg_masked,ConnectDef);

for i = 1:num
    MaxCluster = max(MaxCluster,sum(L(:)==i));
end

function [Y,XYZ] = shiNiftiRead(Img)
V = spm_vol(char(Img));
[Y,XYZ] = spm_read_vols(V);


