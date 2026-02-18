function Gz = shiSpmImgCompress(Img,Gz,do_dt16,do_deleteImg)

% compresses images to .gz
% 2026.1.21

Img = cellfun(char(Img));

if ~exist('Gz','var') || isempty(Gz)
    Gz = [Img{1},'.gz'];
end

if ~exist('dt16','var') || isempty(do_dt16)
    do_dt16 = false;
end

if ~exist('deleteImg','var') || isempty(do_deleteImg)
    do_deleteImg = false;
end

if do_dt16
    xximg16(Img);
end

tar(Gz,Img);

if exist(Gz,'file')
    if do_deleteImg
        cellfun(@delete,Img)
    end
else
    error('cannot find %s',Gz);
end


function xximg16(f)
for i = 1:length(f)
    [pt,~,xt] = fileparts(f{i});
    if ~matches(lower(xt),{'.nii','.img'}), continue; end
    V = spm_vol(f{i});
    if V.dt(1)~=64, continue; end
    Y = single(spm_read_vols(V));
    V.dt(1) = 16;
    tmp = fullfile(pt,['_tmp',sprintf('%+f',randn),xt]);
    movefile(f{i},tmp);
    spm_write_vol(V,Y);
    assert(exist(f{i},'file'),sprintf('fail to convert %s (renamed to %s)',f{i},tmp));
    delete(tmp);
end