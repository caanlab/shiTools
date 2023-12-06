function [X,Roi,RoiFormat] = shiSpmRoiXtr_Nan(Img,varargin)

% extracts summary value(s) of ROI(s) from images, allowing non-existing Img filenames (see shiSpmRoiXtr.m)
%

Img = cellstr(char(Img));
X = nan(size(Img));

idx = cellfun(@(x)exist(x,'file'),Img)>0;

Img1 = Img(idx);
[X1,Roi,RoiFormat] = shiSpmRoiXtr(Img1,varargin{:});

X(idx) = X1;
