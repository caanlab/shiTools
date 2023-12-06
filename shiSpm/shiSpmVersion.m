function Ver = shiSpmVersion(Version)

% returns or switches spm version (8 vs  12)
%
% v = shiSpmVersion        - returns version
% shiSpmVersion('spm12')   - switches to spm12
%
% zhenhao shi

PathSpm = fileparts(which('spm.m'));
[PathRoot,Ver] = fileparts(PathSpm);

if nargin < 1
    return;
end

path_8 = genpath(fullfile(PathRoot,'spm8'));
path_12 = fullfile(PathRoot,'spm12');
path_12_all = genpath(path_12);


if strcmpi(Version,'spm8')
    rmpath(path_12_all);
    addpath(path_8);
elseif strcmpi(Version,'spm12')
    rmpath(path_8);
    addpath(path_12_all);
%     clear classes;
%     addpath(fullfile(path_12,'compat'));
else
    error('input should be either spm8 or spm12');
end

PathSpm = fileparts(which('spm.m'));
[~,Ver] = fileparts(PathSpm);
