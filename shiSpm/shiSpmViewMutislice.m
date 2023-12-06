function datBrain = shiSpmViewMutislice(Underlay,Overlay,VoxThres,SliceCut,Coord,ColorLimit,ColorMap,OverlayInterpolation)

%

% set underlay
if isempty(Underlay) || isequal(lower(Underlay),'default')
    Underlay = which('shi_defaultUnderlay.nii');
else
    Underlay = char(Underlay);
end

% set overlays
Overlay = cellstr(char(Overlay));

% set overlay threshold
if ~exist('VoxThres','var') || isempty(VoxThres)
    VoxThres = -inf(size(Overlay));
elseif isscalar(VoxThres) && isnumeric(VoxThres)
    VoxThres = repmat(VoxThres,length(Overlay),1);
elseif numel(VoxThres) ~= length(Overlay) || ~isnumeric(VoxThres)
    error('VoxThres should be either empty, a scalar, or a vector of the same length as Overlay');
end

% set overlay color limits
no_color = false;
if ~exist('ColorLimit','var') || isempty(ColorLimit)
    no_color = true;
    ColorLimit = [VoxThres(:),ceil(abs(VoxThres(:).*1.5))];
elseif numel(ColorLimit)==2 && isnumeric(ColorLimit)
    ColorLimit = repmat(ColorLimit(:)',length(Overlay),1);
elseif isequal(size(ColorLimit),[length(Overlay),2])
    error('ColorLimit should be either empty, a 1x2 vector, or a Nx2 matrix');
end

% set overlay color map
if ~exist('ColorMap','var') || isempty(ColorMap)
    ColorMap = {autumn,winter,spring,summer,pink,copper,bone};
    ColorMap = ColorMap( mod((1:numel(Overlay))-1,7) + 1 );
end

% check slice cut direction
if ~exist('SliceCut','var') || isempty(SliceCut)
    SliceCut = 'z';
elseif ~ismember(lower(SliceCut),{'x','y','z'})
    error('SliceCut should be either x, y, or z');
end

% check coordinates
if ~exist('Coord','var') || isempty(Coord)
    Coord = -100:5:100;
end

% set interpolation orders
if ~exist('OverlayInterpolation','var') || isempty(OverlayInterpolation)
    OverlayInterpolation = 0;
end

% reslice underlay and overlays to high-resolution
[V,Y] = shiSpmResliceRead(char([{Underlay};Overlay]),OverlayInterpolation);
YUnder = Y(:,:,:,1);
YOver = Y(:,:,:,2:end);
YUnder = ( YUnder - nanmin(YUnder(:)) ) ./ ( nanmax(YUnder(:)) - nanmin(YUnder(:)) );
VoxSize = mean(abs(diag(V(1).mat(1:3,1:3))));

for i = 1:length(VoxThres)
    if VoxThres(i) == -Inf
        VoxThres(i) = min(min(min(YOver(:,:,:,i))));
        if no_color
            ColorLimit(i,:) = [VoxThres(i),ceil(abs(VoxThres(i).*1.5))];
        end
    end
end

% convert MNI coordinates to indices and get actual MNI coordinates
[Index,Coord] = shi_mni2ind(SliceCut,Coord,V(1).mat);
selectIndex = Index>0 & Index<size(YUnder,shiIf( isequal(lower(SliceCut),'x'),1, shiIf( isequal(lower(SliceCut),'y'),2, shiIf( isequal(lower(SliceCut),'z'),3,-1 ) ) ));
Coord = Coord(selectIndex);
Index = Index(selectIndex);

% extract slices from underlay and overlays
switch lower(SliceCut)
    case 'x'
        ViewUnder = squeeze(YUnder(Index,:,:));
        ViewUnder = rot90(permute(ViewUnder,[2,3,1]));
        ViewBorder = [false(size(ViewUnder,1),1,size(ViewUnder,3)),true(size(ViewUnder)),false(size(ViewUnder,1),1,size(ViewUnder,3))];
        ViewUnder = [zeros(size(ViewUnder,1),1,size(ViewUnder,3)),ViewUnder,zeros(size(ViewUnder,1),1,size(ViewUnder,3))];
        ViewOver = cell(length(Overlay),1);
        for i = 1:length(Overlay)
            ViewOver{i} = squeeze(YOver(Index,:,:,i));
            ViewOver{i} = rot90(permute(ViewOver{i},[2,3,1]));
            ViewOver{i}(ViewOver{i}<=VoxThres(i)) = NaN;
            ViewOver{i} = [nan(size(ViewOver{i},1),1,size(ViewOver{i},3)),ViewOver{i},nan(size(ViewOver{i},1),1,size(ViewOver{i},3))];
        end
    case 'y'
        ViewUnder = squeeze(YUnder(:,Index,:));
        ViewUnder = rot90(permute(ViewUnder,[1,3,2]));
        ViewBorder = [false(size(ViewUnder,1),1,size(ViewUnder,3)),true(size(ViewUnder)),false(size(ViewUnder,1),1,size(ViewUnder,3))];
        ViewUnder = [zeros(size(ViewUnder,1),1,size(ViewUnder,3)),ViewUnder,zeros(size(ViewUnder,1),1,size(ViewUnder,3))];
        ViewOver = cell(length(Overlay),1);
        for i = 1:length(Overlay)
            ViewOver{i} = squeeze(YOver(:,Index,:,i));
            ViewOver{i} = rot90(permute(ViewOver{i},[1,3,2]));
            ViewOver{i}(ViewOver{i}<=VoxThres(i)) = NaN;
            ViewOver{i} = [nan(size(ViewOver{i},1),1,size(ViewOver{i},3)),ViewOver{i},nan(size(ViewOver{i},1),1,size(ViewOver{i},3))];
        end
    case 'z'
        ViewUnder = squeeze(YUnder(:,:,Index));
        ViewUnder = fliplr(rot90(ViewUnder));
        ViewBorder = [false(size(ViewUnder,1),1,size(ViewUnder,3)),true(size(ViewUnder)),false(size(ViewUnder,1),1,size(ViewUnder,3))];
        ViewUnder = [zeros(size(ViewUnder,1),1,size(ViewUnder,3)),ViewUnder,zeros(size(ViewUnder,1),1,size(ViewUnder,3))];
        ViewOver = cell(length(Overlay),1);
        for i = 1:length(Overlay)
            ViewOver{i} = squeeze(YOver(:,:,Index,i));
            ViewOver{i}(ViewOver{i}<=VoxThres(i)) = NaN;
            ViewOver{i} = fliplr(rot90(ViewOver{i}));
            ViewOver{i} = [nan(size(ViewOver{i},1),1,size(ViewOver{i},3)),ViewOver{i},nan(size(ViewOver{i},1),1,size(ViewOver{i},3))];
        end
end

% line up slices in one row
ViewUnder = shi_line_up(ViewUnder);
ViewOver = shi_line_up(ViewOver);
ViewBorder = shi_line_up(ViewBorder);

% get overlay transparency and expand to RGB layers
TransOver = cell(size(ViewOver));
for i = 1:numel(TransOver)
    TransOver{i} = repmat(isnan(ViewOver{i}),1,1,3);
end
ViewBorder = repmat(ViewBorder,1,1,3);

% convert value to color
ViewUnder = shi_get_color(ViewUnder,[.2,.8],gray(512)); % set underlay color limits to [.2,.8]
for i = 1:numel(ViewOver)
    ViewOver{i} = shi_get_color(ViewOver{i},ColorLimit(i,:),ColorMap{i});
end

% resize if necessary
if VoxSize>1
    ViewUnder = imresize(ViewUnder,[NaN, round(size(ViewUnder,2)*VoxSize/length(Coord))*length(Coord) ]) + eps;
    ViewBorder = imresize(ViewBorder,[NaN, round(size(ViewBorder,2)*VoxSize/length(Coord))*length(Coord) ],'nearest');
    for i = 1:numel(ViewOver)
        ViewOver{i} = imresize(ViewOver{i},[NaN, round(size(ViewOver{i},2)*VoxSize/length(Coord))*length(Coord) ],'nearest');
        TransOver{i} = imresize(TransOver{i},[NaN, round(size(TransOver{i},2)*VoxSize/length(Coord))*length(Coord) ],'nearest');
    end
end

% merge underlay and overlay
datBrain = ViewUnder;
for i = 1:numel(ViewOver)
    datBrain(~TransOver{i}) = ViewOver{i}(~TransOver{i});
end
datBrain(~ViewBorder) = 0;

% add text
H = 30;
W = size(datBrain,2)/length(Coord);
datText = cell(size(Coord));
for i = 1:length(Coord)
    datText{i} = shi_print_num(Coord(i),H,W);
end
datText = repmat([datText{:}],1,1,3);
datBrain = cat(1,datBrain,datText);

% add color bar
fHi = size(datBrain,1);
if fHi > 60
    datCbar = shi_add_cbar(fHi,ColorLimit,ColorMap);
    datBrain = [datBrain,datCbar];
end

% show
axes('Units','Normalize','Position',[0 0 1 1]);
datBrain = shiPhotoChop_x(datBrain,min(size(datBrain,2),1600),[0 0 0]); % set ouput picture width to 1600 pixels
imshow(datBrain);
figposit = get(gcf,'Position');
set(gcf,'Position',[100,100,1200,figposit(4)./figposit(3).*1200]);

if ~nargout
    datBrain = [];
end

function cbar = shi_add_cbar(h,clim,cmap)
txtHeight = 40;
txtWidth = 150;
barWidth = 30;
flim = cell(size(clim));
for i = 1:size(flim,1)
    flim{i,1} = shi_print_num(clim(i,1),txtHeight,txtWidth,'%.2f',3,8,18);
    flim{i,2} = shi_print_num(clim(i,2),txtHeight,txtWidth,'%.2f',3,18,8);
    flim{i,1} = repmat(flim{i,1},1,1,3) + eps;
    flim{i,2} = repmat(flim{i,2},1,1,3) + eps;
end
h = h - txtHeight - txtHeight;
fmap = cell(size(cmap));
for i = 1:numel(cmap)
    fmap{i} = flip(imresize(permute(cmap{i},[1,3,2]),[h,barWidth],'nearest'));
    ww = txtWidth - barWidth;
    bg1 = zeros(h,round(ww/2),3);
    bg2 = zeros(h,ww-round(ww/2),3);
    fmap{i} = [bg1,fmap{i},bg2];
end
for i = 1:numel(cmap)
    fmap{i} = [flim{i,2};fmap{i};flim{i,1}];
end
cbar = [fmap{:}];



function [Ind,actCoord] = shi_mni2ind(SliceCut,Coord,Mat)
Ind = nan(numel(Coord),1);
actCoord = nan(numel(Coord),1);
for i = 1:length(Coord)
    switch lower(SliceCut)
        case 'x'
            xInd = round([Coord(i),0,0,1]/Mat');
            xCoord = xInd*Mat';
            Ind(i) = xInd(1);
            actCoord(i) = xCoord(1);
        case 'y'
            xInd = round([0,Coord(i),0,1]/Mat');
            xCoord = xInd*Mat';
            Ind(i) = xInd(2);
            actCoord(i) = xCoord(2);
        case 'z'
            xInd = round([0,0,Coord(i),1]/Mat');
            xCoord = xInd*Mat';
            Ind(i) = xInd(3);
            actCoord(i) = xCoord(3);
    end
end



function vv = shi_line_up(v)
if iscell(v)
    for i = 1:numel(v)
        v{i} = shi_line_up(v{i});
    end
    vv = v;
    return;
end
vv = [];
for i = 1:size(v,3)
    vv = [vv,v(:,:,i)]; %#ok<AGROW>
end



function RGB = shi_get_color(Value,ColorLimit,ColorMap)

if ischar(ColorMap)
    ColorMap = eval(ColorMap);
end

Lower = ColorLimit(1);
Upper = ColorLimit(2);
if Upper < Lower
    Lower = ColorLimit(2);
    Upper = ColorLimit(1);
    ColorMap = ColorMap(end:-1:1,:);
end 

Value(Value>=Upper) = Upper;
Value(Value<=Lower) = Lower;
Value(isnan(Value)) = Lower;

ValueColor = round( (Value-Lower) ./ (Upper-Lower) .* (size(ColorMap,1)-1) + 1);
RGB = nan([size(Value),3]);
size(RGB);
for i = 1:size(Value,1)
    for j = 1:size(Value,2)
        RGB(i,j,:) = ColorMap(ValueColor(i,j),:);
    end
end



function ViewNum = shi_print_num(Number,Height,Width,NumFormat,Coat,Hat,Shoe)

if ~exist('NumFormat','var')
    NumFormat = '%+.5g';
end
if ~exist('Coat','var')
    Coat = 3;
end
if ~exist('Hat','var')
    Hat = 3;
end
if ~exist('Shoe','var')
    Shoe = 3;
end
Char = sprintf(NumFormat,Number);

Singleton = cell(size(Char));
for i = 1:length(Char)
    Singleton{i} = get_singleton(Char(i));
    Singleton{i} = add_coat(Singleton{i},Coat);
    Singleton{i} = add_hat(Singleton{i},Hat);
    Singleton{i} = add_shoe(Singleton{i},Shoe);
end

ViewNum = [Singleton{:}];
ViewNum = change_size(ViewNum,Height,Width);



function vv = change_size(v,hh,ww)
[h,w] = size(v);
if h/w >= hh/ww
    vv = make_fatter(v,hh,ww);
else
    vv = make_taller(v,hh,ww);
end



function vv = make_fatter(v,hh,ww)
vv = imresize(v,[hh,NaN]);
ww = ww - size(vv,2);
bg1 = false(hh,round(ww/2));
bg2 = false(hh,ww-round(ww/2));
vv = [bg1,vv,bg2];



function vv = make_taller(v,hh,ww)
vv = imresize(v,[NaN,ww]);
hh = hh - size(vv,1);
bg1 = false(round(hh/2),ww);
bg2 = false(hh-round(hh/2),ww);
vv = [bg1;vv;bg2];



function Singleton = add_coat(Singleton,Coat)
Singleton = [zeros(size(Singleton,1),Coat),Singleton,zeros(size(Singleton,1),Coat)];



function Singleton = add_hat(Singleton,Hat)
Singleton = [zeros(Hat,size(Singleton,2));Singleton];



function Singleton = add_shoe(Singleton,Shoe)
Singleton = [Singleton;zeros(Shoe,size(Singleton,2))];



function Singleton = get_singleton(Digit)
sz = [39,20];
x{1} = [1:39:742,2:39:743,3:39:744,4:39:745,5:39:746];
x{2} = [1:22,40:61,79:100,118:139,157:178];
x{3} = x{2}+585;
x{4} = x{1}+17;
x{5} = x{2}+17;
x{6} = x{3}+17;
x{7} = x{4}+17;
x{8} = [96:39:681,97:39:682,98:39:683,99:39:684]-78;
x{9} = [324:339,363:378,402:417,441:456]-78;
x{10} = [307:312,346:351,385:390,424:429,463:468,502:507];

switch Digit
    case '0'
        Ind = setdiff(1:7,4);
    case '1'
        Ind = [3,6];
    case '2'
        Ind = setdiff(1:7,[2,6]);
    case '3'
        Ind = setdiff(1:7,[2,5]);
    case '4'
        Ind = [2,3,4,6];
    case '5'
        Ind = setdiff(1:7,[3,5]);
    case '6'
        Ind = setdiff(1:7,3);
    case '7'
        Ind = [1,3,6];
    case '8'
        Ind = 1:7;
    case '9'
        Ind = setdiff(1:7,5);
    case '+'
        Ind = 8:9;
    case '-'
        Ind = 8;
    case '.'
        Ind = 10;
end

Singleton = false(sz);
Singleton([x{Ind}]) = true;

