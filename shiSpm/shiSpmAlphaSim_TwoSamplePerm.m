function [MaxClusterSize_Pos,MaxClusterSize_Neg,MaxT_Pos,MaxT_Neg] = shiSpmAlphaSim_TwoSamplePerm(Img1,Img2, WorkingDir, MaskImg, OutTxt, ConnectDef, Voxel_P, n_Iteration)

% performs permutation-based Monte-Carlo simulation using the "Eklund" approach
%
% (masking procedure will be revised)


figure(173);
pause(0.01);

Img1 = cellstr(char(Img1));
Img2 = cellstr(char(Img2));
ImgAll = [Img1;Img2];

%% reslice mask
PermDir = shiMkdir(fullfile(WorkingDir,'PermDir'));
MaskImg_Resliced = char(shiSpmReslice(MaskImg,ImgAll{1},WorkingDir,0,'perm_RelicedMask_'));
Y_MaskImg_Resliced = shiNiftiRead(MaskImg_Resliced) > 0;

%% threshold p to t
DF = length(ImgAll)-2;
Voxel_T = tinv(1-Voxel_P,DF);

%% flip coins
Cov = nan(size(ImgAll));
CovRaw = [ones(size(Img1));2*ones(size(Img2))];
for iter = 1:n_Iteration
    Cov(:,iter) = CovRaw(randperm(length(ImgAll)));
end

%% permute
MaxClusterSize_Pos = nan(n_Iteration,numel(Voxel_P));
MaxClusterSize_Neg = nan(n_Iteration,numel(Voxel_P));
MaxT_Pos = nan(n_Iteration,1);
MaxT_Neg = nan(n_Iteration,1);

for iter = 1:n_Iteration
    
    shiDisp({'Two Sample T Permutation:',sprintf('%d/%d',iter,n_Iteration)});

    subplot(9,2,[1,2]);
    pcolor([repmat(Cov(:,iter)',[2,1]),[NaN;NaN]]);
    axis off

    shiSpmAlphaSim_TwoSamplePerm_TwoSampleT(PermDir,ImgAll(Cov(:,iter)==1),ImgAll(Cov(:,iter)==2));
    if exist(fullfile(PermDir,'spmT_0001.nii'),'file')
        Y_TImg = shiNiftiRead(fullfile(PermDir,'spmT_0001.nii'));
    else
        Y_TImg = shiNiftiRead(fullfile(PermDir,'spmT_0001.img'));
    end
    Y_TImg(~Y_MaskImg_Resliced) = 0;

    for p = 1:numel(Voxel_T)
        Y_TImg_masked_Pos = ( Y_TImg > Voxel_T(p) ) + 0;
        Y_TImg_masked_Neg = ( Y_TImg < -Voxel_T(p) ) + 0;
        MaxClusterSize_Pos(iter,p) = shiSpmAlphaSim_TwoSamplePerm_MaxCluster(Y_TImg_masked_Pos,ConnectDef);
        MaxClusterSize_Neg(iter,p) = shiSpmAlphaSim_TwoSamplePerm_MaxCluster(Y_TImg_masked_Neg,ConnectDef);
    end
    MaxT_Pos(iter) = nanmax(Y_TImg(:));
    MaxT_Neg(iter) = -nanmin(Y_TImg(:));
    
    subplot(9,2,3:2:9);
    MaxClusterSize_Pos_sorted = sort(MaxClusterSize_Pos(1:iter,:),1,'descend');
    semilogx(MaxClusterSize_Pos_sorted);
    hold on;
    semilogx([0.05*iter,0.05*iter],[0,nanmax(MaxClusterSize_Pos(:))],'k');
    if iter > 1
        legend(shiStrConcat('p<',num2str( Voxel_P(:) ,'%-g'),', k_{0.05}Pos=',num2str(ceil(eps+interp1(MaxClusterSize_Pos_sorted,0.05*iter,'linear','extrap'))','%-d')));
    end
    hold off;

    subplot(9,2,11:2:17);
    MaxT_Pos_sorted = sort(MaxT_Pos(1:iter,:),1,'descend');
    semilogx(MaxT_Pos_sorted);
    hold on;
    semilogx([0.05*iter,0.05*iter],[0,nanmax(MaxT_Pos(:))],'k');
    if iter > 1
        legend(['T_{0.05}Pos=',num2str(interp1(MaxT_Pos_sorted,0.05*iter,'linear','extrap'),'%-.2f')]);
    end
    hold off;

    subplot(9,2,4:2:10);
    MaxClusterSize_Neg_sorted = sort(MaxClusterSize_Neg(1:iter,:),1,'descend');
    semilogx(MaxClusterSize_Neg_sorted);
    hold on;
    semilogx([0.05*iter,0.05*iter],[0,nanmax(MaxClusterSize_Neg(:))],'k');
    if iter > 1
        legend(shiStrConcat('p<',num2str( Voxel_P(:) ,'%-g'),', k_{0.05}Neg=',num2str(ceil(eps+interp1(MaxClusterSize_Neg_sorted,0.05*iter,'linear','extrap'))','%-d')));
    end
    hold off;

    subplot(9,2,12:2:18);
    MaxT_Neg_sorted = sort(MaxT_Neg(1:iter,:),1,'descend');
    semilogx(MaxT_Neg_sorted);
    hold on;
    semilogx([0.05*iter,0.05*iter],[0,nanmax(MaxT_Neg(:))],'k');
    if iter > 1
        legend(['T_{0.05}Neg=',num2str(interp1(MaxT_Neg_sorted,0.05*iter,'linear','extrap'),'%-.2f')]);
    end
    hold off;

    shiSpmAlphaSim_TwoSamplePerm_ResetWorkingDir(PermDir);
    
    pause(0.001);

end

cd(WorkingDir);

%%
MaxClusterSize_sorted_Pos = MaxClusterSize_Pos;
MaxClusterSize_sorted_Neg = MaxClusterSize_Neg;
fid = fopen(OutTxt,'a+');
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
dlmwrite(OutTxt,MaxClusterSize_Pos_sorted,'-append','delimiter','\t');




function shiSpmAlphaSim_TwoSamplePerm_ResetWorkingDir(PermDir)
%%
delete(fullfile(char(PermDir),'*'));



function shiSpmAlphaSim_TwoSamplePerm_TwoSampleT(PermDir,Img_1,Img_2)
%%

matlabbatch{1}.spm.stats.factorial_design.dir = {PermDir};
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = Img_1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = Img_2;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
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
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '01_G1_G2';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

spm_jobman('serial',matlabbatch);



function MaxCluster = shiSpmAlphaSim_TwoSamplePerm_MaxCluster(Y_TImg_masked,ConnectDef)
%%
MaxCluster = 0;

[L,num] = spm_bwlabel(Y_TImg_masked,ConnectDef);

for i = 1:num
    MaxCluster = max(MaxCluster,sum(L(:)==i));
end

function [Y,XYZ] = shiNiftiRead(Img)
V = spm_vol(char(Img));
[Y,XYZ] = spm_read_vols(V);
