function [MaxClusterSize,MaxT] = shiSpmAlphaSim_OneSamplePerm_par(Img, WorkingDir, MaskImg, OutTxt, ConnectDef, Voxel_P, n_Iteration)

% performs permutation-based Monte-Carlo simulation using the "Eklund" approach (parfor version)
%

Img = cellstr(char(Img));
MaskImg = char(shiFullFileName(MaskImg));

%% reslice mask
PermDir_prefix = fullfile(WorkingDir,'PermDir');
% MaskImg_Resliced = char(shiSpmReslice(MaskImg,Img{1},WorkingDir,0,'perm_RelicedMask_'));
% Y_MaskImg_Resliced = shiNiftiRead(MaskImg_Resliced) > 0;

[~,Y_MaskImg_Resliced,~,~,~] = shiSpmMaskRead(Img{1},MaskImg);

%% threshold p to t
DF = length(Img)-1;
Voxel_T = tinv(1-Voxel_P,DF);

%% flip coins
Cov = nan(length(Img),n_Iteration);
PermDir = cell(n_Iteration,1);
for iter = 1:n_Iteration
    Cov(:,iter) = 2*unidrnd(2,length(Img),1)-3;
    PermDir{iter,1} = [PermDir_prefix,'_',num2str(iter),'_',num2str(unidrnd(1000000))];
end

%% permute
MaxClusterSize = nan(n_Iteration,numel(Voxel_P));
MaxT = nan(n_Iteration,1);

parfor iter = 1:n_Iteration
    
%     shiDisp({'One Sample T Permutation:',sprintf('%d/%d',iter,n_Iteration)});
%     
%     subplot(911);
%     pcolor([repmat(Cov(:,iter)',[2,1]),[NaN;NaN]]);
%     axis off

    shiSpmAlphaSim_OneSamplePerm_OneSampleRegress(PermDir{iter},Img,Cov(:,iter));
    if exist(fullfile(PermDir{iter},'spmT_0001.nii'),'file')
        Y_TImg = shiNiftiRead(fullfile(PermDir{iter},'spmT_0001.nii'));
    else
        Y_TImg = shiNiftiRead(fullfile(PermDir{iter},'spmT_0001.img'));
    end
    Y_TImg(~Y_MaskImg_Resliced) = 0;

    tmp_K = [];
    for p = 1:numel(Voxel_T)
        Y_TImg_masked = ( Y_TImg > Voxel_T(p) ) + 0;
        tmp_K(1,p) = shiSpmAlphaSim_OneSamplePerm_MaxCluster(Y_TImg_masked,ConnectDef);
    end
    MaxClusterSize(iter,:) = tmp_K;
    MaxT(iter) = nanmax(Y_TImg(:));
    
%     subplot(9,1,2:5);
%     MaxClusterSize_sorted = sort(MaxClusterSize(1:iter,:),1,'descend');
%     semilogx(MaxClusterSize_sorted);
%     hold on;
%     semilogx([0.05*iter,0.05*iter],[0,nanmax(MaxClusterSize(:))],'k');
%     if iter > 1
%         legend(shiStrConcat('p<',num2str( Voxel_P(:) ,'%-g'),', k_{0.05}=',num2str(ceil(eps+interp1(MaxClusterSize_sorted,0.05*iter,'linear','extrap'))','%-d')));
%     end
%     hold off;
% 
%     subplot(9,1,6:9);
%     MaxT_sorted = sort(MaxT(1:iter,:),1,'descend');
%     semilogx(MaxT_sorted);
%     hold on;
%     semilogx([0.05*iter,0.05*iter],[0,nanmax(MaxT(:))],'k');
%     if iter > 1
%         legend(['T_{0.05}=',num2str(interp1(MaxT_sorted,0.05*iter,'linear','extrap'),'%-.2f')]);
%     end
%     hold off;

    shiSpmAlphaSim_OneSamplePerm_ResetWorkingDir(PermDir{iter});
    
%     pause(0.001);

end

cd(WorkingDir);

%%
MaxClusterSize_sorted = MaxClusterSize;
fid = fopen(OutTxt,'a+');
fprintf(fid,'\n\n%s\n%s\n\n',repmat('#',1,19),shiTime(0));
for p = 1:numel(Voxel_T)
    MaxClusterSize_sorted(:,p) = sort(MaxClusterSize_sorted(:,p),'descend');
    fprintf(fid,'p_voxel<%g\t',Voxel_P(p));
end
fprintf(fid,'Tmax\tp_Cluster\n\n');
fclose(fid);
MaxClusterSize_sorted = [MaxClusterSize_sorted';sort(MaxT,'descend')';(1:n_Iteration)/n_Iteration]';
dlmwrite(OutTxt,MaxClusterSize_sorted,'-append','delimiter','\t');




function shiSpmAlphaSim_OneSamplePerm_ResetWorkingDir(PermDir)
%%
delete(fullfile(char(PermDir),'*'));




function shiSpmAlphaSim_OneSamplePerm_OneSampleRegress(PermDir,Img,Cov)
%%
matlabbatch{1}.spm.stats.factorial_design.dir = {char(PermDir)};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = Img;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
matlabbatch{1}.spm.stats.factorial_design.cov.c = Cov;
matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'OneSampleReg';
matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'c';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('serial',matlabbatch);



function MaxCluster = shiSpmAlphaSim_OneSamplePerm_MaxCluster(Y_TImg_masked,ConnectDef)
%%
MaxCluster = 0;

[L,num] = spm_bwlabel(Y_TImg_masked,ConnectDef);

for i = 1:num
    MaxCluster = max(MaxCluster,sum(L(:)==i));
end

function [Y,XYZ] = shiNiftiRead(Img)
V = spm_vol(char(Img));
[Y,XYZ] = spm_read_vols(V);

