function outImg = shiSpm3dTo4d(Img3d,Img4d)
 
% combines 3D nifti images into a single 4D image
%
% outImg = shiSpm3dTo4d(Img3d,Img4d)
% 
%   Img3d - a cell array of 3D nifti file names
%   Img4d - a string of output 4D file name
%   outImg - cell of output 4D image file name (with full path)
% 
%    ###########
% by Zhenhao Shi @ 2018-9-11
%    ###########
% 

[Path4d,~,~] = fileparts(Img4d);

if isempty(Path4d)
    Img4d = fullfile(pwd,Img4d);
elseif ~isempty(Path4d) && ~isfolder(Path4d)
    mkdir(Path4d);
end

matlabbatch{1}.spm.util.cat.vols = Img3d;
matlabbatch{1}.spm.util.cat.name = Img4d;
matlabbatch{1}.spm.util.cat.dtype = 0;

spm_jobman('serial',matlabbatch);
outImg = shiFullFileName(Img4d);

