function shiSpmImgCalc(InputImg,FuncHandle,OutputImg,Mask,dt)

% performs any voxel-wise calculations on nifti images
%
% shiSpmImgCalc(InputImg,FuncHandle,OutputImg)
% shiSpmImgCalc(InputImg,FuncHandle,OutputImg,Mask)
% shiSpmImgCalc(InputImg,FuncHandle,OutputImg,Mask,dt)
%   can do complicated calculations on nifti images. Theoretically, this
%   function can run all types of vector computations that MATLAB is
%   capable of.
% 
%   InputImg   - a string or a cell array of strings for input image file
%                name(s)
%   FuncHandle - handle for a vector-in vector-out function. the input
%                vector and output vector should be of equal length as
%                InputImg and OutputImg, respectively
%   OutputImg  - a string or a cell array of strings for output image file
%                name(s)
%   Mask       - a string of inclusive mask file, or '' (default) to 
%                not use mask
%   dt         - output format (default = 64) (see spm_type)
% 
% Example:
%    InputImg  = {
%                'A1.img';
%                'A2.img';
%                'A3.img';
%                'A4.img';
%                'A5.img';
%                'A6.img';
%                'A7.img';
%                'B1.img';
%                'B2.img';
%                'B3.img';
%                'B4.img';
%                'B5.img';
%                'B6.img';
%                'B7.img';
%                };
%    OutputImg = {
%                'Result/Spearman_R.img'; % spearman corr between A and B
%                'Result/Pearson_R.img'; % pearson corr between A and B
%                };
%    FuncHandle = @(x)([corr(x(1:7),x(8:14),'type','spearman'),corr(x(1:7),x(8:14))]);
%    shiSpmImgCalc(InputImg,FuncHandle,OutputImg);
% 
%    ###########
% by Zhenhao Shi @ 2018-04-11
%    ###########
% 

if ~exist('Mask','var') || isempty(Mask)
    applyMask = false;
else
    applyMask = true;
end

if ~exist('dt','var') || isempty(dt)
    try
        DT = spm_vol(char(InputImg));
        dt = cell(length(DT),1);
        [dt{:}] = deal(DT.dt);
        dt=cell2mat(dt);
        dt = max(dt(:,1));
    catch
        warning('cannot read "dt" info');
        dt = 64;
    end
end

InputImg = char(InputImg);
OutputImg = cellstr(char(OutputImg));

[OutputImgPath,~,~] = shiFileParts(OutputImg);
for i = 1:length(OutputImgPath)
    if ~isempty(OutputImgPath{i}) && ~isfolder(OutputImgPath{i})
        mkdir(OutputImgPath{i});
    end
end


[V_InputImg,Y_InputImg,~] = shiSpmResliceRead(InputImg);
if applyMask
    [~,Y_Mask] = shiSpmMaskRead(deblank(InputImg(1,:)),Mask,'incl');
else
    Y_Mask = true(size(Y_InputImg(:,:,:,1)));
end


Y_OutputImg = nan([size(Y_InputImg(:,:,:,1)),length(OutputImg)]);

V_OutputImg = cell(length(OutputImg),1);
for o = 1:length(V_OutputImg)
    V_OutputImg{o} = struct('fname',   OutputImg{o},...
                            'dim',     V_InputImg(1).dim(1:3),...
                            'dt',      [dt spm_platform('bigend')],...
                            'mat',     V_InputImg(1).mat,...
                            'descrip', ['shiSpmImgCalc:',func2str(FuncHandle)]);
end

fprintf('shiSpmImgCalc: slice   0');
for k = 1:size(Y_InputImg,3)
    fprintf('\b\b\b%3d',k);
    for j = 1:size(Y_InputImg,2)
        for i = 1:size(Y_InputImg,1)
            
            if ~Y_Mask(i,j,k)
                continue;
            end
            
            Input = squeeze(Y_InputImg(i,j,k,:));
            Output = FuncHandle(Input);
            Y_OutputImg(i,j,k,:) = Output;

        end
    end
end




for o = 1:length(V_OutputImg)
    spm_write_vol(V_OutputImg{o},Y_OutputImg(:,:,:,o));
end
