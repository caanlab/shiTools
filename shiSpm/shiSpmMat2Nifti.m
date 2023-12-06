function [fname,V] = shiSpmMat2Nifti(matrix,fname)

% writes numeric matrices as nifti images
%
% Zhenhao Shi 2020-4-16
%

if ~iscell(matrix)
    matrix = {matrix};
end

fname = cellstr(char(fname));

if numel(matrix) ~= numel(fname)
    error('number of matrices must be equal to number of output nii filenames');
end

for i = 1:numel(matrix)
    if ndims(matrix{i})>3 || ndims(matrix{i})<1
        error('matrix must be 1-d, 2-d or 3-d');
    end
end

V = cell(numel(matrix),1);

for i = 1:numel(matrix)
    V{i} = spm_vol;
    V{i}(1).fname = fname{i};
    V{i}.pinfo = [1;0];
    V{i}.mat = eye(4);
    [V{i}.dim(1),V{i}.dim(2),V{i}.dim(3)] = size(matrix{i});
    spm_write_vol(V{i},matrix{i});
end
