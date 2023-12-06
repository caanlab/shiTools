function Nii = shiSpmImg2Nii(Img,doDeleteImg)

% converts .img/hdr file pairs to nii files
% 
% Nii = shiSpmImg2Nii(Img)
% 
%    ###########
% by Zhenhao Shi @ 2015-3-26
%    ###########

Img = cellstr(char(Img));

if ~exist('doDeleteImg','var') || isempty(doDeleteImg)
    doDeleteImg = false;
end

Nii = cell(size(Img));

for j = 1:length(Img)
    fprintf('% 5d %s',j,Img{j});
    [path,name,ext] = fileparts(Img{j});
    if ~strcmpi(ext,'.img')
        warning('  is not a .img file\n');
        continue;
    end
    Nii{j} = fullfile(path,[name,'.nii']);
    V = spm_vol(Img{j});
    Y = spm_read_vols(V);
    if ndims(Y) == 4
        warning('  is 4-D and cannot be converted\n');
        continue;
    end
    V.fname = Nii{j};
    V = rmfield(V,'private');
    spm_write_vol(V,Y);
    fprintf(' --> %s\n',Nii{j});
    if doDeleteImg
        delete(Img{j});
        delete([Img{j}(1:end-4),'.hdr']);
    end
end