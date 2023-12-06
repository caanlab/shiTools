function shiSpmJob(SpmBatch)

% executes saved SPM batch file or struct
%
% shiSpmJob(SpmBatch)
%
%   SpmBatch - string of cell array for batch file name, or job struct
%
% requires
%   SPM 8/12
% 
%    ###########
% by Zhenhao Shi @ 2021-4-5
%    ###########
% 

if isstruct(SpmBatch) || (iscell(SpmBatch) && all(cellfun(@isstruct,SpmBatch)))
    spm_jobman('serial',SpmBatch);
else
    SpmBatch = cellstr(char(SpmBatch));
    for i = 1:length(SpmBatch)
        spm_jobman('serial',SpmBatch{i});
    end
end