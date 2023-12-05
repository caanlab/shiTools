function MatOut = shiMatShift(MatIn,Background,Direction)

% moves matrix along directions, with size unchanged, showing background near borders
%
% MatOut = shiMatShift(MatIn,Background,Direction)
%   Direction - [moveUpDown,moveLeftRight,movePrevioussliceNextslice,...]
%
% zhenhao shi

MatOut = MatIn;
MatOut(:) = Background;

indIn = cell(ndims(MatIn),1);
indOut = cell(ndims(MatIn),1);

for i = 1:ndims(MatIn)
    indIn{i} = 1:size(MatIn,i);
    indOut{i} = 1:size(MatIn,i);
end
    
for i = 1:numel(Direction)
%     max(1,1-Direction(i))
%     min(size(MatIn,i),size(MatIn,i)-Direction(i))
%     max(1,1+Direction(i))
%     min(size(MatIn,i),size(MatIn,i)+Direction(i))
    indIn{i} = max(1,1-Direction(i)):min(size(MatIn,i),size(MatIn,i)-Direction(i));
    indOut{i} = max(1,1+Direction(i)):min(size(MatIn,i),size(MatIn,i)+Direction(i));
end

MatOut(indOut{:}) = MatIn(indIn{:});
    