function shiSound(style)

% plays one of six sounds
%
% shiSound(style)
%   style = 1 or 2 or 3 or 4 or 5 or 6
%
%    ###########
% by Zhenhao Shi @ 2014-12-23
%    ###########

if nargin<1 || ~isnumeric(style) || style<1 || style>6 || style~=int32(style)
    style = 1;
end

AudioDir = fullfile(toolboxdir('matlab'),'audiovideo');

switch style
    case 1
        load(fullfile(AudioDir,'gong.mat'));
        sound(y);
    case 2
        load(fullfile(AudioDir,'chirp.mat'));
        sound(y);
    case 3
        load(fullfile(AudioDir,'train.mat'));
        sound(y);
    case 4
        load(fullfile(AudioDir,'splat.mat'));
        sound(y);
    case 5
        load(fullfile(AudioDir,'handel.mat'));
        sound(y,Fs*1.35);
    case 6
        load(fullfile(AudioDir,'laughter.mat'));
        sound(y,Fs*1.35);
end

        
