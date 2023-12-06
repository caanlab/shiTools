function outImg = shiSpm4dTo3d(Img4d,Path3d)

% converts 4D nifti image to 3D images
% 
% outImg = shiSpm4dTo3d(Img4d)
% outImg = shiSpm4dTo3d(Img4d,Path3d)
% 
%   Img4d  - a string of input 4D file name. output 3D files will be saved
%            in the same directory of the 4D file
%   outImg - cell of output 4D image file name (with full path)
%   Path3d - where 3D images are to be saved (default: same as Img4d)
% 
%    ###########
% by Zhenhao Shi @ 2015-1-5
%    ###########
% 

if ~exist(Img4d,'file')
    error('cannot find %s',Img4d);
end

Vd = spm_vol(Img4d);
[Path,Name,Ext]=fileparts(Img4d);
if isempty(Path)
    Path=pwd;
end

if nargin < 2
    Path3d = Path;
end
shiMkdir(Path3d);

outImg = cell(length(Vd),1);

for i=1:length(Vd)
    
    W = Vd(i);
    X = spm_read_vols(Vd(i));
    
    imgnum = num2str(W.n(1));
    numb = '00000';
    numb(end-length(imgnum)+1:end) = imgnum;
    
%     W = rmfield(W,'private');
    W.n = [1 1];
    W.fname = fullfile(Path3d,[Name,'_',numb,Ext]);
    outImg{i,1} = W.fname;
    spm_write_vol(W,X);

end
