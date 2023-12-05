function wb_command(varargin)

PathWb = fullfile(fileparts(shiTools),'workbench',shiIf(ismac,'bin_macosx64',shiIf(ispc,'bin_windows64','error')));
WB = fullfile(PathWb, 'wb_command');

STR = '';
for i =  1:nargin
    if iscell(varargin{i})
        for k = 1:length(varargin{i})
            STR = [STR,' ',varargin{i}{k}];
        end
    elseif ischar(varargin{i})
        STR = [STR,' ',varargin{i}];
    elseif isnumeric(varargin{i}) && isscalar(varargin{i})
        STR = [STR,' ',num2str(varargin{i})];
    else
        error('wb_command: %d-th argument not recognized',i);
    end
end

eval(['!',WB,STR]);