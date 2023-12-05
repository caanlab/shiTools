function Photo_out = shiPhotoChop_x(Photo,width,bgcolor,Photo_out)

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
    if ~isequal(im(:,i,:),repmat(bgcolor,size(im,1),1,1))
        markContent = [markContent,i]; %#ok<AGROW>
    end
end

markStart = [markContent([true,diff(markContent)~=1]),size(im,2)+1];

row = 0;
lab = nan(1,size(im,2));
upper = 1;

while any(isnan(lab))
    row = row + 1;
    upper_bef = upper;
    upper = max(markStart(markStart<=width+upper));
    if upper == upper_bef
        error('choose larger width');
    end
    lab(sum(~isnan(lab))+1:upper-1) = row;
end

ROW = cell(row,1);

for i = 1:row
    ROW{i} = im(:,lab==i,:);
    ROW{i}(:,sum(lab==i)+1:width,:) = repmat(bgcolor,size(im,1),width-sum(lab==i),1);
end

oum = cat(1,ROW{:});

if ischar(Photo) && ( ~exist('Photo_out','var') || isempty(Photo_out) )
    [pth,nme,ext] = fileparts(Photo);
    Photo_out = fullfile(pth,[nme,'_chopped',ext]);
end
if ischar(Photo)
    imwrite(oum,Photo_out);
else
    Photo_out = oum;
end