function Photo_out = shiPhotoOverlay(Underlay,Overlay,bgcolor,Photo_out)

under = imread(Underlay);
over = imread(Overlay);

if ~isequal(size(under),size(over))
    error('unequal image sizes');
end

if ~exist('bgcolor','var') || isempty(bgcolor)
    bgcolor = over(1,1,:);
else
    bgcolor = reshape(bgcolor(:),1,1,3);
end

for i = 1:size(over,1)
    for j = 1:size(over,2)
        if isequal(over(i,j,:),bgcolor)
            over(i,j,:) = NaN;
        end
    end
end

under(~isnan(over)) = over(~isnan(over));

if ~exist('Photo_out','var') || isempty(Photo_out)
    [pth,nme,ext] = fileparts(Underlay);
    [~,nme2] = fileparts(Overlay);
    Photo_out = fullfile(pth,[nme,'_overlaid_by_',nme2,ext]);
end
imwrite(under,Photo_out);