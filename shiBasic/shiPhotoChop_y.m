function Photo_out = shiPhotoChop_y(Photo,height,bgcolor,Photo_out)

if ischar(Photo)
    im = imread(Photo);
else
    im = Photo;
end

if ~exist('bgcolor','var') || isempty(bgcolor)
    bgcolor = im(1,1,:);
else
    bgcolor = reshape(bgcolor(:),1,1,3);
end

markContent = [];
for i = 1:size(im,2)
    if ~isequal(im(i,:,:),repmat(bgcolor,1,size(im,2),1))
        markContent = [markContent;i]; %#ok<AGROW>
    end
end

markStart = [markContent([true;diff(markContent)~=1]);size(im,1)+1];

col = 0;
lab = nan(size(im,1),1);
upper = 1;

while any(isnan(lab))
    col = col + 1;
    upper_bef = upper;
    upper = max(markStart(markStart<=height+upper));
    if upper == upper_bef
        error('choose larger width');
    end
    lab(sum(~isnan(lab))+1:upper-1) = col;
end

COL = cell(1,col);

for i = 1:col
    COL{i} = im(lab==i,:,:);
    COL{i}(sum(lab==i)+1:height,:,:) = repmat(bgcolor,height-sum(lab==i),size(im,2),1);
end

oum = cat(2,COL{:});

if ischar(Photo) && ( ~exist('Photo_out','var') || isempty(Photo_out) )
    [pth,nme,ext] = fileparts(Photo);
    Photo_out = fullfile(pth,[nme,'_chopped',ext]);
end
if ischar(Photo)
    imwrite(oum,Photo_out);
else
    Photo_out = oum;
end