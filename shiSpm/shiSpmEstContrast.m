function [SpmFile,Contrast,ContrastName,ContrastType,matlabbatch] = shiSpmEstContrast(SpmFile,Contrast,ContrastName,ContrastType,DeleteExistingContrast,AddZeroBefore,AddZeroAfter)

% reads, adds or rewrites contrast in SPM.mat
%
% [SpmFile,Contrast,ContrastName,ContrastType] = shiSpmEstContrast(SpmFile)
% [SpmFile,Contrast,ContrastName,ContrastType] = shiSpmEstContrast(SpmFile,Contrast)
% [SpmFile,Contrast,ContrastName,ContrastType] = shiSpmEstContrast(SpmFile,Contrast,ContrastName)
% [SpmFile,Contrast,ContrastName,ContrastType] = shiSpmEstContrast(SpmFile,Contrast,ContrastName,ContrastType)
% [SpmFile,Contrast,ContrastName,ContrastType] = shiSpmEstContrast(SpmFile,Contrast,ContrastName,ContrastType,DeleteExistingContrast)
% [SpmFile,Contrast,ContrastName,ContrastType] = shiSpmEstContrast(SpmFile,Contrast,ContrastName,ContrastType,DeleteExistingContrast,AddZeroBefore
% [SpmFile,Contrast,ContrastName,ContrastType] = shiSpmEstContrast(SpmFile,Contrast,ContrastName,ContrastType,DeleteExistingContrast,AddZeroBefore,AddZeroAfter)
%
%  SpmFile                   - SPM.mat file name
%  Contrast                  - cell or matrix (t-con only) of contrasts
%  ContrastName              - contrast names
%  ContrastType              - 'T' or 'F' contrasts (Contrast must be cell for 'F' contrasts)
%  DeleteExistingContrast    - default: 0
%  AddZeroBefore             - add ? number of zeros before each element of Contrast (e.g. in the case of parametric modulation)
%  AddZeroAfter              - add ? number of zeros after each element of Contrast (e.g. in the case of derivatives)
%
% Zhenhao Shi 2023/10/04
%

if ~exist('Contrast','var') || isempty(Contrast)
    try
        SPM = load(SpmFile);
        Contrast = cellfun(@transpose,{SPM.SPM.xCon.c}','UniformOutput',false);
        ContrastName = {SPM.SPM.xCon.name}';
        ContrastType = {SPM.xCon.STAT}';
    catch
        [Contrast,ContrastName,ContrastType] = deal({});
    end
    return;
end

if ~exist('DeleteExistingContrast','var') || isempty(DeleteExistingContrast)
    DeleteExistingContrast = 0;
end

if ~exist('AddZeroBefore','var') || isempty(AddZeroBefore)
    AddZeroBefore = 0;
end
if ~exist('AddZeroAfter','var') || isempty(AddZeroAfter)
    AddZeroAfter = 0;
end

if isnumeric(Contrast)
    tmpContrast = Contrast;
    Contrast = cell(size(Contrast,1),1);
    for i = 1:numel(Contrast)
        Contrast{i} = tmpContrast(i,:);
    end
end
Contrast = Contrast(:);

if ~exist('ContrastName','var') || isempty(ContrastName)
    ContrastName = cell(numel(Contrast),1);
    for i = 1:numel(Contrast)
        ContrastName{i} = sprintf('Contrast%02d',i);
    end
end
ContrastName = cellstr(char(ContrastName));

if length(Contrast) ~= length(ContrastName)
    error('contrast name length must equal contrast length');
end

if AddZeroBefore+AddZeroAfter > 0
    
    tmpContrast = Contrast;
    for i = 1:length(Contrast)
        Contrast{i} = zeros(1,(1+AddZeroBefore+AddZeroAfter)*numel(tmpContrast{i}));
        ind = (1+AddZeroBefore):(1+AddZeroBefore+AddZeroAfter):(1+AddZeroBefore+AddZeroAfter)*numel(tmpContrast{i});
        Contrast{i}(ind) = tmpContrast{i}(:);
    end
    
end
        
if ~exist('ContrastType','var') || isempty(ContrastType) || (ischar(ContrastType) && strcmpi(ContrastType,'T'))
    ContrastType = repmat({'T'},length(Contrast),1);
elseif ischar(ContrastType) && strcmpi(ContrastType,'F')
    ContrastType = repmat({'F'},length(Contrast),1);
end
ContrastType = cellstr(char(ContrastType));

%-----------------------------------------------------------------------
% Job saved on 16-Dec-2019 15:54:20 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

matlabbatch{1}.spm.stats.con.spmmat = {SpmFile};
for i = 1:numel(Contrast)
    if strcmpi(ContrastType{i},'T')
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = ContrastName{i};
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.weights = Contrast{i};
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'replsc';
    elseif strcmpi(ContrastType{i},'F')
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.name = ContrastName{i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.weights = Contrast{i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.sessrep = 'replsc';
    else
        error('unrecognized contrast type #%d: %s',i,ContrastType{i});
    end
end
matlabbatch{1}.spm.stats.con.delete = DeleteExistingContrast*1;

spm_jobman('serial',matlabbatch);

if nargout < 2
    return;
end

SPM = load(SpmFile);
Contrast = cellfun(@transpose,{SPM.SPM.xCon.c}','UniformOutput',false);
ContrastName = {SPM.SPM.xCon.name}';
ContrastType = {SPM.SPM.xCon.STAT}';
