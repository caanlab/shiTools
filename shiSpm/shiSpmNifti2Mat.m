function [matrix,V,XYZ] = shiSpmNifti2Mat(fname)

% reads nifti images into numeric matrices
%
% Zhenhao Shi 2020-4-16
%

fname = cellstr(char(fname));

matrix = cell(numel(fname),1);
V = cell(size(matrix));
XYZ = cell(size(matrix));

for i = 1:numel(matrix)
    V{i} = spm_vol(fname{i});
    [matrix{i},XYZ{i}] = spm_read_vols(V{i});
    if ndims(matrix{i})==4
        warning('%s is a 4-d image',fname{i});
    end
end
