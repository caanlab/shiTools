function matlabbatch = shiSpmStatRegress(Dir,Img,Vector,VectorName,doContrast,Mask,write4D)

% conducts whole-brain regression analysis
% 
% shiSpmStatRegress(Dir,Img,Vector)
% shiSpmStatRegress(Dir,Img,Vector,VectorName)
%   Dir        - string, where results are to be saved
%   Img        - input nifti images
%   Vector     - a matrix of the same number of rows as Img to regress
%                voxelwisely against Img. Each column is a variable
%   VectorName - cell array of strings for variable names (default =
%                'Vector1', 'Vector2', ...)
%   doContrast - whether to do [1] and [-1] contrasts for each predictor
%                (default = true)
% 
%    ###########
% by Zhenhao Shi @ 2024-6-24
%    ###########
% 

PWD = pwd;

shiMkdir(Dir);
Img=cellstr(char(Img));

if size(Img,1) ~= size(Vector,1)
    error('unmatched observation number');
end

Img_orig = Img;
Vector_orig = Vector; %#ok<*NASGU>
[Img,Vector,AnyMiss] = shi_deNaN(Img,Vector);

if ~exist('VectorName','var') || isempty(VectorName)
    for j = 1:size(Vector,2)
        VectorName{j} = ['Vector',num2str(j)];
    end
else
    VectorName = cellstr(char(VectorName));
    if size(Vector,2) ~= numel(VectorName)
        error('wrong number of names for vectors');
    end
end

if size(Img,1) ~= size(Vector,1)
    error('unmatched observation number');
end

if ~exist('Mask','var') || isempty(Mask)
    Mask = {''};
else
    Mask = cellstr(char(Mask));
    if numel(Mask)~=1
        error('Mask should be either left empty or specified as one inclusive mask image filename');
    end
end


if AnyMiss
    warning('\n\n  #############################\n  ##                         ##\n  ##       - WARNING -       ##\n  ##                         ##\n  ##  Missing values found!  ##\n  ##                         ##\n  #############################\n  ##  %s\n\n',Dir);
    save(fullfile(Dir,'_MissingValuePresent.mat'),'Img_orig','Img','Vector_orig','Vector');
end


matlabbatch{1}.spm.stats.factorial_design.dir = {Dir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = Img;


for j = 1:size(Vector,2)
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).c = Vector(:,j);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).cname = VectorName{j};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j).iCC = 1;
end
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).c = ones(length(Img),1);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).cname = 'Constant';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(j+1).iCC = 5;




matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = Mask;
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

if ~exist('doContrast','var') || isempty(doContrast) || doContrast

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
    matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
    matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');

    for j = 1:size(Vector,2)
        matlabbatch{3}.spm.stats.con.consess{j*2-1}.tcon.name = [num2str(j*2-1,'%.2d'),'_Pos_',VectorName{j}];
        matlabbatch{3}.spm.stats.con.consess{j*2-1}.tcon.convec = [zeros(1,j-1),1,zeros(1,size(Vector,2)-j),0];
        matlabbatch{3}.spm.stats.con.consess{j*2-1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{j*2}.tcon.name = [num2str(j*2,'%.2d'),'_Neg_',VectorName{j}];
        matlabbatch{3}.spm.stats.con.consess{j*2}.tcon.convec = [zeros(1,j-1),-1,zeros(1,size(Vector,2)-j),0];
        matlabbatch{3}.spm.stats.con.consess{j*2}.tcon.sessrep = 'none';
    end
    matlabbatch{3}.spm.stats.con.consess{j*2+1}.tcon.name = [num2str(j*2+1,'%.2d'),'_Pos_Constant'];
    matlabbatch{3}.spm.stats.con.consess{j*2+1}.tcon.convec = [zeros(1,j),1];
    matlabbatch{3}.spm.stats.con.consess{j*2+1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{j*2+2}.tcon.name = [num2str(j*2+2,'%.2d'),'_Neg_Constant'];
    matlabbatch{3}.spm.stats.con.consess{j*2+2}.tcon.convec = [zeros(1,j),-1];
    matlabbatch{3}.spm.stats.con.consess{j*2+2}.tcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.delete = 0;

end

spm_jobman('serial',matlabbatch);


%%
StatType = 'Regress';
Time = shiTime;
if ~exist('write4D','var') || isempty(write4D) || ~write4D
else
    shiSpm3dTo4d(Img,fullfile(Dir,'ConImg4d.img'));
end
% Img = 'ConImg4d.img';
if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','VectorName','Mask','StatType','PWD','Time','matlabbatch','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Dir','Img','Vector','VectorName','Mask','StatType','PWD','Time','matlabbatch');
end


