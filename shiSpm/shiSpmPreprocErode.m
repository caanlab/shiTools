function outImg = shiSpmPreprocErode(TisDepImg,Keep,existAction)

% erodes tissue probability images
%
% outImg = shiSpmPreprocErode(TisDepImg,Keep)
% outImg = shiSpmPreprocErode(TisDepImg,Keep,existAction)
%
% Img        - tissue depth map(s) (see shiSpmPreprocTissueDepth.m)   
% Keep       - percentage of top-prob tissue to be kept (default Keep=0.1, to keep top 10%)
% outImg     - output image filenames: c1_TisDepImg, c2_TisDepImg, c3_TisDepImg
%
% Zhenhao Shi 2019/12/10
%

TisDepImg = cellstr(char(TisDepImg));

if ~exist('Keep','var') || isempty(Keep)
    Keep = .1;
end
if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end

if length(TisDepImg)>1
    outImg = cell(size(TisDepImg));
    for i = 1:length(TisDepImg)
        [outImg{i,1}] = shiSpmPreprocErode(TisDepImg{i},Keep,existAction);
    end
    return;
end

TisDepImg = char(TisDepImg);
V = spm_vol(TisDepImg);
[pth,nme,ext] = fileparts(V.fname);
outImg = {
    fullfile(pth,['ec1_',nme,ext]);
    fullfile(pth,['ec2_',nme,ext]);
    fullfile(pth,['ec3_',nme,ext]);
    };

if any(cellfun(@(x)exist(x,'file'),outImg)) && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('one or more of %s, %s and %s already exists. overwrite?  [y/n] \nK>> ',outImg{3},outImg{2},outImg{3}), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end


%%

Y = spm_read_vols(V);
[e1,e2,e3] = deal(zeros(size(Y)));

Yc1 = Y(Y(:)>=1 & Y(:)<2);
Yc2 = Y(Y(:)>=2 & Y(:)<3);
Yc3 = Y(Y(:)>=3 & Y(:)<4);
Ycut1 = prctile(Yc1,(1-Keep)*100);
Ycut2 = prctile(Yc2,(1-Keep)*100);
Ycut3 = prctile(Yc3,(1-Keep)*100);

if Ycut1 == max(Yc1)
    e1(Y>=Ycut1 & Y<2) = 1;
else
    e1(Y>Ycut1 & Y<2) = 1;
end
if Ycut2 == max(Yc2)
    e2(Y>=Ycut2 & Y<3) = 1;
else
    e2(Y>Ycut2 & Y<3) = 1;
end
if Ycut3 == max(Yc3)
    e3(Y>=Ycut3 & Y<4) = 1;
else
    e3(Y>Ycut3 & Y<4) = 1;
end


V.fname = outImg{1};
V.descrip = [V.descrip, sprintf('; eroded c1, Type=%s', Keep)];
spm_write_vol(V,e1);

V.fname = outImg{2};
V.descrip = [V.descrip, sprintf('; eroded c2, Type=%s', Keep)];
spm_write_vol(V,e2);

V.fname = outImg{3};
V.descrip = [V.descrip, sprintf('; eroded c3, Type=%s', Keep)];
spm_write_vol(V,e3);


