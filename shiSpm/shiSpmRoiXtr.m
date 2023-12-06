function [X,Roi,RoiFormat] = shiSpmRoiXtr(Img,Roi,Radius,SummFunc)

% extracts summary value(s) of ROI(s) from images
%
% Roi can be a binary .img/.nii image, or MNI coordinates ([x,y,z]) or MNI coordinates+radius ([x,y,z,radius]), or a structure (see shiSpmRoiFormat)
% The variable "Radius" is used only when no radius is given in "Roi"
% 
% shiSpmRoiXtr(Img,Roi)
% shiSpmRoiXtr(Img,Roi,Radius)
% shiSpmRoiXtr(Img,Roi,SummFunc)
% shiSpmRoiXtr(Img,Roi,Radius,SummFunc)% X: n_Img-by-n_Roi matrix
%
% Example:
% X=shiSpmRoiXtr({'con_0001.nii';'con_0002.nii'},{which('shi_AAL_RIns.img'),[0 50 5 3],[45 2 2 5]},'mean'}
%   returns the mean voxel values at the AAL-defined right insula, the
%   5-mm-radius spheric region at 0/50/5, and the 3-mm-radius spheric
%   region at 45/2/2, from the images con_0001.nii and con_0002.nii.
%
%    ###########
% by Zhenhao Shi @ 2021-3-22
%    ###########
% 

if isempty(Roi)
    X = [];
    Roi = [];
    RoiFormat = [];
    return;
end

thres = 0;

if nargin == 3 && ischar(Radius)
    SummFunc = Radius; % shiSpmRoiXtr(Img,Roi,SummFunc)
end

if ~exist('Radius','var') || isempty(Radius) || ischar(Radius) || ~(Radius>=0)
    Radius = 5;
end


%% check Img 

File_Img = char(Img);
for i = size(File_Img,1)
    if ~exist(deblank(File_Img(i,:)),'file')
        error('cannot find %s',File_Img(i,:))
    end
end

[Roi,RoiFormat] = shiSpmRoiFormat(Roi,Radius);


%%

V_Img = spm_vol(File_Img);
X = NaN(size(File_Img,1),length(Roi));

for r = 1:length(Roi)

    switch RoiFormat{r}
        
        case 'Mask'

            File_Roi = Roi{r};
            V_Roi = spm_vol(File_Roi);
            [Y_Roi,XYZ_Roi] = spm_read_vols(V_Roi);
            XYZ_Roi = XYZ_Roi(:, Y_Roi > thres);
            
        case 'Atlas'
            
            File_Roi = Roi{r}.atlas;
            V_Roi = spm_vol(File_Roi);
            [Y_Roi,XYZ_Roi] = spm_read_vols(V_Roi);
            XYZ_Roi = XYZ_Roi(:, Y_Roi == Roi{r}.value);
            
        case 'Sphere'
            
            XYZ_Roi = shi_read_img_xyz(V_Img(1));
            XYZ_Roi = XYZ_Roi( : , sqrt( sum( ( XYZ_Roi' - repmat(Roi{r}(1:3),size(XYZ_Roi,2),1) ).^2 , 2 ) ) <= Roi{r}(4) );
            
        otherwise
            
            error('Roi # %d: unknown Roi format',r);
            
    end
    
    XYZ_Roi_voxel = (V_Img(1).mat) \ [XYZ_Roi' ones(size(XYZ_Roi,2),1)]';
    xY_Img = spm_get_data(File_Img,XYZ_Roi_voxel);

    if ~exist('SummFunc','var') || isempty(SummFunc) || strcmpi(SummFunc(1:3),'mea')
        X(:,r) = nanmean(xY_Img,2);
        
    elseif strcmpi(SummFunc(1:3),'med')
        X(:,r) = nanmedian(xY_Img,2);
        
    elseif strcmpi(SummFunc(1:3),'max')
        X(:,r) = nanmax(xY_Img,[],2);
        
    elseif strcmpi(SummFunc(1:3),'min')
        X(:,r) = nanmin(xY_Img,[],2);
        
    elseif strcmpi(SummFunc(1:3),'eig') % see spm_regions
        try
            not_nan = all(~isnan(xY_Img));
            xY_Img = xY_Img(:,not_nan);
            [m,n]   = size(xY_Img);
            if m > n
                [~,s,v] = svd(xY_Img'*xY_Img);
                s       = diag(s);
                v       = v(:,1);
                u       = xY_Img*v/sqrt(s(1));
            else
                [~,s,u] = svd(xY_Img*xY_Img');
                s       = diag(s);
                u       = u(:,1);
                v       = xY_Img'*u/sqrt(s(1));
            end
            d       = sign(sum(v));
            u       = u*d;
            X(:,r)       = u*sqrt(s(1)/n);
        catch
            X(:,r) = NaN;
        end
        
    else
        warning('unknown summarizing methods. use ''mean'' instead');
        X(:,r) = nanmean(xY_Img,2);
        
    end

end


function XYZ_Img = shi_read_img_xyz(V_Img)
[R,C,P] = ndgrid(1:V_Img(1).dim(1), 1:V_Img(1).dim(2), 1:V_Img(1).dim(3));
RCP     = [R(:)'; C(:)'; P(:)'; ones(1,numel(R))];
XYZ_Img   = V_Img(1).mat(1:3,:) * RCP;