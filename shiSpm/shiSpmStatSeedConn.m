function shiSpmStatSeedConn(Dir,Img,Roi,Mask,write4D)
 
% conduts seed-based connectivity analysis.
% 
% shiSpmStatSeedConn(Dir,Img,Roi)
%   Dir         - string, path to save the results
%   Img         - cell array of strings, .img/.nii data file names
%   Roi         - string of ROI file name (.img/.nii), or MNI coordinates
%                 or coordinates+radius (see shiSpmRoiXtr)
% 
%    ###########
% by Zhenhao Shi @ 2020-6-29
%    ###########
%

X=shiSpmRoiXtr(Img,Roi);

if ~exist('write4D','var') || isempty(write4D)
    write4D = false;
end

if ~exist('Mask','var') || isempty(Mask)
    Mask = {''};
else
    Mask = cellstr(char(Mask));
    if numel(Mask)~=1
        error('Mask should be either left empty or specified as one inclusive mask image filename');
    end
end

shiSpmStatCorr(Dir,Img,X,Mask,write4D);


if exist(fullfile(Dir,'StatInfo.mat'),'file')
    save(fullfile(Dir,'StatInfo.mat'),'Roi','-append');
else
    save(fullfile(Dir,'StatInfo.mat'),'Roi');
end

