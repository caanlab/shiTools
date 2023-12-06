function shiSpmResetOrigin(Images)

% resets origin of nifti images [0,0,0] to the center of the space
%
% thanks to Kanchana Jagannathan

Images = cellstr(char(Images));

for i=1:numel(Images)
    fprintf('%s\n',Images{i});
    V    = spm_vol(Images{i});
    M    = V.mat;
    vox  = sqrt(sum(M(1:3,1:3).^2));
    if det(M(1:3,1:3))<0
        vox(1) = -vox(1);
    end
    orig = (V.dim(1:3)+1)/2;
    off  = -vox.*orig;
    M    = [vox(1) 0      0      off(1)
            0      vox(2) 0      off(2)
            0      0      vox(3) off(3)
            0      0      0      1     ];
    spm_get_space(Images{i},M);
end

