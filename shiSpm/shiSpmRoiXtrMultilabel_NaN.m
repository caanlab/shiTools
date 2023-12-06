function [X,NodeValue,NodeName] = shiSpmRoiXtrMultilabel_NaN(Img,varargin)

% extracts summary value(s) of ROIs defined by an atlas, allowing non-existing Img filenames (see shiSpmRoiXtrMultilabel.m)
%

Img = cellstr(char(Img));

idx = cellfun(@(x)exist(x,'file'),Img)>0;

Img1 = Img(idx);
[X1,NodeValue,NodeName] = shiSpmRoiXtrMultilabel(Img1,varargin{:});

X = nan(size(Img,1),size(X1,2));

X(idx,:) = X1;
