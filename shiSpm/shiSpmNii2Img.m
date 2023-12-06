function Img = shiSpmNii2Img(Nii)

% converts .nii files to .img/.hdr file pairs
% 
% Nii = shiSpmNii2Img(Img)
% 
%    ###########
% by Zhenhao Shi @ 2015-3-26
%    ###########

Nii = cellstr(char(Nii));

cnt = 0;
Img = {};

for i = 1:length(Nii)
    nii = shiFullFileName(Nii{i});
    for j = 1:length(nii)
        fprintf('%s\n',nii{j});
        [path,name,ext] = fileparts(nii{j});
        if ~strcmpi(ext,'.nii')
            warning('%s is not a .nii file\n',nii{j});
            continue;
        end
        cnt = cnt+1;
        Img{cnt,1} = fullfile(path,[name,'.img']);
        V = spm_vol(nii{j});
        Y = spm_read_vols(V);
        V.fname = [name,'.img'];
        V.private.dat.fname = V.fname;
        cwd = pwd;
        cd(path);
        spm_write_vol(V,Y);
        cd(cwd);
    end;
end;