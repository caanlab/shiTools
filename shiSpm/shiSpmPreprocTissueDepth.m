function [DepthImg,LabelImg] = shiSpmPreprocTissueDepth(c1,c2,c3,existAction)

% uses tissue probability maps (c1,c2,c3) to calculate multilabel tissue atlas and tissue depth map (e.g. depth value of 1.23 indicates the voxel being buried in tissue c1 deeper than 23% of other c1 voxels)
%
% Zhenhao Shi 2020-7-23

c1 = cellstr(char(c1));
c2 = cellstr(char(c2));
c3 = cellstr(char(c3));

if length(c1) ~= length(c2) || length(c1) ~= length(c3)
    error('input tissue images must be of the same length');
end

if length(c1) > 1
    DepthImg = cell(size(c1));
    LabelImg = cell(size(c1));
    for i = 1:length(c1)
        [DepthImg{i,1},LabelImg{i,1}] = shiSpmPreprocTissueDepth(c1{i},c2{i},c3{i});
    end
    return;
end

[pth, nme, ext] = fileparts(char(c1));
V = spm_vol(char(c1));

V_dep = V;
V_lab = V;
V_dep.fname = fullfile(pth,['TisDep_',nme,ext]);
V_lab.fname = fullfile(pth,['TisLab_',nme,ext]);
V_dep.dt(1) = 64;
V_lab.dt(1) = 64;
DepthImg = V_dep.fname;
LabelImg = V_lab.fname;

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end
if exist(DepthImg,'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',DepthImg), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end


%%

[Y,XYZ] = spm_read_vols(spm_vol(char([c1;c2;c3])));
Y(:,:,:,4) = 1 - sum(Y,4);
[~,Label] = max(Y,[],4); % (Label: 1=c1, 2=c2, 3=c3, 4=other)
Label(Label==4) = 0;

msk = cell(3,1);
border = cell(3,1);
Depth = zeros(size(Label));
for c = 1:3
    msk{c} = Label==c;
    border{c} = sBorder(msk{c});
    msk{c} = find(msk{c});
    chunksize = round(sqrt(length(msk{c})));
    nchunk = ceil(length(msk{c})/chunksize);
    fprintf('\ncomputing tissue depth c%d:   0%%',c);
    dist = [];
    for ch = 1:nchunk
        fprintf('\b\b\b\b%3d%%',floor(ch/nchunk*100));
        beg = (ch-1)*chunksize+1;
        fin = min(length(msk{c}),ch*chunksize);
        dist(beg:fin,1) = min( pdist2(XYZ(:,msk{c}(beg:fin))',XYZ(:,border{c}(:))') , [] , 2 );
    end
    dist = (tiedrank(dist)-0.5)/length(dist);
    Depth(msk{c}) = c+dist;
end

fprintf('\n');
spm_write_vol(V_dep,Depth);
spm_write_vol(V_lab,Label);


function border = sBorder(msk)

Shift = [
    +1 00 00
    -1 00 00
    00 +1 00
    00 -1 00
    00 00 +1
    00 00 -1
    ];

msk_shift = false([size(msk),size(Shift,1)]);

for i = 1:size(Shift,1)
    msk_shift(:,:,:,i) = shiMatShift(msk,false,Shift(i,:));
end

border = any(msk_shift,4) & ~msk;
