function ImgMasked = shiSpmMask(ImgRaw,Mask,Type,Interpolation,Prefix,SaveMask)

% performs inclusive or exclusive masking
%
% ImgMasked = shiSpmMask(ImgRaw,Mask)
% ImgMasked = shiSpmMask(ImgRaw,Mask,Type)
% ImgMasked = shiSpmMask(ImgRaw,Mask,Type,Interpolation)
% ImgMasked = shiSpmMask(ImgRaw,Mask,Type,Interpolation,SaveMask)
%   returns masked images
%   all masked images are to be resliced to the first raw image
%   all masks are to be resliced to the first mask before calculation
%
%   ImgRaw          - raw unmasked image file(s), or value [1] to calculate mask only
%   Mask            - mask image file(s)
%   Type            - cellstr of two strings. First string indicates 'inclusive' (default) or 'exclusive'. Second string indicates 'AND' (default) or 'OR logic
%                   - note that "inclusive"/"exclusive" is computed AFTER "AND"/"OR", see below:
%                           if {'inclusive','AND'}, then mask is calculated as:     Mask_1 & Mask_2 & Mask_3 ...
%                           if {'inclusive','OR'},  then mask is calculated as:     Mask_1 | Mask_2 | Mask_3 ...
%                           if {'exclusive','AND'}, then mask is calculated as:   ~(Mask_1 & Mask_2 & Mask_3 ...)
%                           if {'exclusive','OR'},  then mask is calculated as:   ~(Mask_1 | Mask_2 | Mask_3 ...)
%   Prefix          - prefix for output masked image (default = 'Masked_')
%   SaveMask        - if true, the mask will be saved as ['SaveMask_',ImgRaw{i}]
% see spm_mask.m

if ~exist('SaveMask','var') || isempty(SaveMask)
    SaveMask = false;
end

if ~exist('Type','var') || isempty(Type)
    Type = {'inclusive','AND'};
elseif ischar(Type)
    if strcmpi(Type(1:3),'inc') || strcmpi(Type(1:3),'exc')
        Type = {Type,'AND'};
    elseif strcmpi(Type,'AND') || strcmpi(Type,'OR')
        Type = {'inclusive',Type};
    end
elseif iscellstr(Type) && length(Type) == 1
    Type = char(Type);
    if strcmpi(Type(1:3),'inc') || strcmpi(Type(1:3),'exc')
        Type = {Type,'AND'};
    elseif strcmpi(Type,'AND') || strcmpi(Type,'OR')
        Type = {'inclusive',Type};
    else
        error('''Type'' must indicate inclusive/exclusive and AND/OR (e.g. {''incl'',''AND''');
    end
elseif (~strcmpi(Type{1}(1:3),'inc') && ~strcmpi(Type{1}(1:3),'exc')) || (~strcmpi(Type{2},'AND') && ~strcmpi(Type{2},'OR'))
    error('''Type'' must indicate inclusive/exclusive and AND/OR (e.g. {''incl'',''AND''');
end

if ~exist('Interpolation','var') || isempty(Interpolation)
    Interpolation = 0;
end

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'Masked_';
end

P_Mask = char(Mask);
if isequal(ImgRaw,1)
    ImgRaw = deblank(P_Mask(1,:));
    SaveMaskOnly = true;
    SaveMask = true;
else
    SaveMaskOnly = false;
end
P_Img = char(ImgRaw);

if SaveMaskOnly
    ImgMasked = '';
else
    if ischar(ImgRaw) && size(ImgRaw,1) == 1
        [pt,nme,xt] = fileparts(ImgRaw);
        ImgMasked = fullfile(pt,[Prefix,nme,xt]);
    else
        [pt,nme,xt] = shiFileParts(cellstr(char(ImgRaw)));
        ImgMasked = shiStrConcat(pt,filesep,Prefix,nme,xt);
    end
end

% if ~nargin
%     [P_Mask, sts] = spm_select(Inf,'image','Images to compute mask from');
%     if ~sts, return; end
%     [P_Img, sts] = spm_select(Inf,'image','Images to apply mask to');
%     if ~sts, return; end
% end

V_Mask = spm_vol(P_Mask);
V_Img = spm_vol(P_Img);

thresh = zeros(numel(V_Mask),1);

m_Mask = numel(V_Mask);
m_Img = numel(V_Img);

%-Create headers
%--------------------------------------------------------------------------
V_Out = V_Img;
for i=1:m_Img
    [pth,nm,ext,num] = spm_fileparts(deblank(V_Out(i).fname));
    V_Out(i).fname      = fullfile(pth,[Prefix, nm, ext, num]);
    V_Out(i).descrip    = 'Masked';
    V_Out(i).mat        = V_Out(1).mat;
    V_Out(i).dim(1:3)   = V_Out(1).dim(1:3);
    if isfield(V_Img(i),'descrip')
        V_Out(i).descrip    = [V_Img(i).descrip,'; Masked'];
    end
end
if ~SaveMaskOnly
    V_Out  = spm_create_vol(V_Out);
end
M   = V_Out(1).mat;
dim = V_Out(1).dim(1:3);

%-Compute masked images
%--------------------------------------------------------------------------

Y_SaveMask = false(dim);

for j=1:dim(3)

    if strcmpi(Type{2},'AND')
        msk = true(dim(1:2));
    elseif strcmpi(Type{2},'OR')
        msk = false(dim(1:2));
    end
    Mi  = spm_matrix([0 0 j]);

    % Load slice j from all images
    for i=1:m_Mask
        M1  = M\V_Mask(i).mat\Mi;
        %if sum((M1(:)-Mi(:)).^2<eps) M1 = Mi; end;

        img = spm_slice_vol(V_Mask(i),M1,dim(1:2),[0 NaN]);
        msk = msk & isfinite(img);

        if ~spm_type(V_Mask(i).dt(1),'nanrep')
            if strcmpi(Type{2},'AND')
                msk = msk & (img ~= 0);
            elseif strcmpi(Type{2},'OR')
                msk = msk | (img ~= 0);
            end
        end

        if strcmpi(Type{2},'AND')
            msk = msk & (img > thresh(i));
        elseif strcmpi(Type{2},'OR')
            msk = msk | (img > thresh(i));
        end

    end

    if strcmpi(Type{1}(1:3),'exc')
        msk = ~msk;
    end

    % Write the images
    if ~SaveMaskOnly
        for i=1:m_Img
            M1        = M\V_Img(i).mat\Mi;
            img       = spm_slice_vol(V_Img(i),M1,dim(1:2),Interpolation);
            img(~msk) = NaN;
            V_Out(i)     = spm_write_plane(V_Out(i),img,j);
        end
    end

    Y_SaveMask(:,:,j) = msk;

end

if SaveMask
    V_SaveMask = V_Out(1);
    [pth,nm,ext,num] = spm_fileparts(V_SaveMask.fname);
    V_SaveMask.fname = fullfile(pth,['SaveMask_', nm, ext, num]);
    V_SaveMask.dt = [2 0];
    V_SaveMask.descrip = 'SaveMask';
    spm_write_vol(V_SaveMask,Y_SaveMask);
end

