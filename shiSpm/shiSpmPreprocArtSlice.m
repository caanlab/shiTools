function [outImg,ArtSliceLogTxt,ArtSliceSummTxt] = shiSpmPreprocArtSlice(Img,Thres,Prefix,existAction)

% runs art_slice function of ARTRepair

% The bad slices show up best in the raw images, and the bad slice data 
% spreads around in slicetiming and realignment. So check for bad slices 
% before any SPM preprocessing. It is suggested to use this function only 
% with very noisy data.

Img = cellstr(char(Img));
[Pth,Nme,Ext] = shiFileParts(Img);
outImg = shiStrConcat(Pth,filesep,Prefix,Nme,Ext);

if ~exist('existAction','var') || isempty(existAction)
    existAction = 'ask';
elseif ~ismember(lower(existAction),{'ask','overwrite'})
    error('input not recognized');
end
if exist(outImg{1},'file') && strcmpi(existAction,'ask')
    ACTION = input(strrep(sprintf('%s already exists. overwrite?  [y/n] \nK>> ',outImg{1}), '\', '\\'),'s');
    switch lower(ACTION)
        case 'n'
            error('aborted');
        case 'y'
            warning('overwritting...');
        otherwise
            error('input not recognized');
    end
end


if ~exist('Thres','var') || isempty(Thres)
    Thres = 18;
    %  Threshold above sample means to filter slices (def=18)
    %  15 is very visible on contrast image. 8 is slightly visible.
end

repair_flag = 1; % default: Repair Bad Slices and Write BadSliceLog ( repairs only the bad slices )
mask_flag = 1; % default: Automatic ( will generate ArtifactMask image )

if ~exist('Prefix','var') || isempty(Prefix)
    Prefix = 'g'; % default='g' if repair_flag and mask_flag are default
end

shi_art_slice(spm_vol(char(Img)),Thres,repair_flag,mask_flag,Prefix);


ArtSliceLogTxt = fullfile(Pth,['artslice_log_',Nme{1},'.txt']);
ArtSliceSummTxt = fullfile(Pth,['artslice_summ_',Nme{1},'.txt']);


function shi_art_slice(Pimages, OUTSLICE, repair_flag, mask_flag,   Prefix)

[~,tmp_img1_name] = fileparts(Pimages(1).fname);

% FORMAT art_slice  (v3.1)
%
% SUMMARY
%   Reads a set of images and writes a new set of images
%   after filtering the data for noise.  User is asked to choose the
%   repair method or methods. This program is best applied
%   to the raw images, so the cleaned output images can be fed into
%   slicetiming or realignment algorithms.
%
%   V.2 uses the scale factor of input images for the output images
%   to improve compatibility when output is used by external programs
%   (SPM and FSL use different default scaling.)
%
%   To automatically screen an image data set for bad slices, use the
%   Detect and Repair Bad Slices option which writes a BadSliceLog
%   of all the bad slices detected and repaired.
%
%   No original data is removed- the old images are all preserved.
%
% User is asked via GUI to supply:
%   - Set of images
%   - A repair method
%   - If repairing bad slices, user can choose a threshold.
% POSSIBLE REPAIRS
% 1. Filter all the images with a 3-point median filter in time. An
%     excellent filter for block designs with TR=2 or less. If TR > 2 sec.,
%     or for Rapid Event
%     experiments where events are separated by 2 TR's or less, this filter
%     may reduce sensitivity and clobber activations.
%     Adds a prefix "f" to the cleaned output images.
% 2. Detect and repair bad slices. Derives new values for the bad slice
%     using linear interpolation of the before and after volumes.
%     If the same slice is bad in multiple volumes, uses spatial
%     interpolation with neighboring slices instead of pure temporal
%     interpolation (added logic for v3).
%     Bad slices are detected when the amount of data scattered outside
%     the head is at least T above the usual amount for the slice. (The usual
%     amount is determined as the average of the best two of the first three
%     volumes. It differs for each slice.) The default
%     value of T is 5, which the user can adjust. A preview estimate
%     of the amount of repaired data at T=5 is shown in the Matlab window.
%     This filter removes outliers, but may reduce activations, so it's
%     safest when fewer than 5% of the slices are cleaned up.
%     Adds a prefix "g" to the cleaned output images.
%     Writes a text file BadSliceLog of all bad slices detected.
% 3. Eliminate data outside the head. Generates a head mask automatically,
%     and writes it as ArtifactMask.img. Sets voxels outside the head
%     to zero. This process may help realignment to be more accurate when
%     there is large clutter in the field of view, e.g. prenatal imaging.
%     Adds a prefix "h" to the cleaned output images.
% 4. Combinations of methods 1 and 3, or 2 and 3.
%
% NOTES
%  1. Use art_movie to compare the input and repaired images.
%  2. First and last volumes will not be filtered.
% Paul Mazaika - December 2004.  V.2  August 2006,
% v2.2  allows user to specify a head mask instead of automatic. July 2007
% v2.3  compatible with multiple spm versions.  Mar09
% v3.0  added code for cases when slice used for interpolation is
%       artifactual. Dorian Pustina, Oct 2012
% v3.1  supports SPM12 Dec2014 pkm.

% Configure while preserving old SPM versions
spmv = spm('Ver'); spm_ver = 'spm5';  % chooses spm_select to read vols
if (strcmp(spmv,'SPM2')) spm_ver = 'spm2'; end
if (strcmp(spmv,'SPM2') || strcmp(spmv,'SPM5')) spm_defaults;
else spm('Defaults','fmri'); end


% PARAMETERS
OUTSLICEdef =  18;  %  Threshold above sample means to filter slices
%  15 is very visible on contrast image. 8 is slightly visible.

% THREEE CONTROL OPTIONS
%   Mask the images    - Set amask = 1.
%   Median filter all the images  - Set allfilt = 1.
%   Detect and repair bad slices  - Set slicefilt = 1.
if strcmp(spm_ver,'spm5')
    %mask = spm_select(1, 'image', 'Select mask image in functional space');
    if exist('Pimages','var') == 0, Pimages  = spm_select(Inf,'image','select images'); end
else  % spm2 version
    %mask = spm_get(1, '.img', 'Select mask image in functional space');
    if exist('Pimages','var') == 0, Pimages  = spm_get(Inf,'image','select images'); end;
end

P = spm_vol(Pimages);

if exist('repair_flag','var') == 0
    repair_flag = spm_input('Which repair methods to use?', 1, 'm', ...
        ' Repair Bad Slices and Write BadSliceLog ( repairs only the bad slices ) | Median Filter All Data  ( aggressive filter for very noisy data ) | Eliminate data outside head | Repair Bad Slices and Eliminate data outside head | Median Filter All Data and Eliminate data outside head',...
        [1 2 3 4 5], 1);
end

if ( repair_flag == 3 ) amask = 1; allfilt = 0; slicefilt = 0; end
if ( repair_flag == 1 ) amask = 0; allfilt = 0; slicefilt = 1; end
if ( repair_flag == 2 ) amask = 0; allfilt = 1; slicefilt = 0; end
if ( repair_flag == 4 ) amask = 1; allfilt = 0; slicefilt = 1; end
if ( repair_flag == 5 ) amask = 1; allfilt = 1; slicefilt = 0; end

%  Set the prefix correctly on the written files
% % % % % % % if  ( amask ==1 & allfilt == 0 & slicefilt == 0 ) prechar = 'h'; end
% % % % % % % if  ( amask ==1 & allfilt == 1 & slicefilt == 0 ) prechar = 'fh'; end
% % % % % % % if  ( amask ==1 & allfilt == 0 & slicefilt == 1 ) prechar = 'gh'; end
% % % % % % % if  ( amask ==0 & allfilt == 1 & slicefilt == 0 ) prechar = 'f'; end
% % % % % % % if  ( amask ==0 & allfilt == 0 & slicefilt == 1 ) prechar = 'g'; end
prechar = Prefix;

fprintf('\n NEW IMAGE FILES WILL BE CREATED');
fprintf('\n The filtered scan data will be saved in the same directory');
fprintf('\n with %s pre-pended to their filenames.\n',prechar);

if amask == 1 | slicefilt == 1 %  Automask options
    if exist('mask_flag','var') == 0
        mask_flag = spm_input('Which mask image to use?', 1, 'm', ...
            'Automatic ( will generate ArtifactMask image ) | User specified mask ',...
            [1 2], 1);
    end
    if mask_flag == 2
        if strcmp(spm_ver,'spm5')
            maskimg = spm_select(1, '.img', 'Select mask image in functional space');
        else  % spm2 version
            maskimg = spm_get(1, '.img', 'Select mask image in functional space');
        end
        Automask = spm_read_vols(spm_vol(maskimg));
        maskcount = sum(sum(sum(Automask)));  %  Number of voxels in mask.
        voxelcount = prod(size(Automask));    %  Number of voxels in 3D volume.
    else  %  mask_flag == 1
        fprintf('\n Generated mask image is written.\n');
        fprintf('\n');
        Pnames = P(1).fname;
        Automask = shi_art_automask(Pnames(1,:),-1,1);
        maskcount = sum(sum(sum(Automask)));  %  Number of voxels in mask.
        voxelcount = prod(size(Automask));    %  Number of voxels in 3D volume.
    end
end

if amask == 1
    fprintf('\n ELIMINATE DATA OUTSIDE HEAD');
    fprintf('\n All scans will have voxels outside the head set to zero.\n');
end

if slicefilt == 1  % Prepare some thresholds for slice testing.
    fprintf('\n INTERPOLATE THROUGH BAD SLICES\n');
    % Find the slice orientation used to collect the data
    [ vx, vy, vz ] = size(Automask);
    orient = 0;
    if ( vx < vy & vx < vz ) orient = 1; disp(' Remove bad Sagittal slices'); end
    if ( vy < vx & vy < vz ) orient = 2; disp(' Remove bad Coronal slices'); end
    if ( vz < vx & vz < vy ) orient = 3; disp(' Remove bad Axial slices'); end
    nslice = min([vx vy vz]);
    if ( orient == 0 )
        disp('Cannot determine slice orientation for bad slice filtering.')
        return;
    end
    % Find 3 samples of slice baseline activity outside the head.
    p = zeros(3,nslice);
    for i = 1:3
        Y1 = spm_read_vols(P(i));
        Y1 = ( 1 - Automask ).*Y1;
        % Get the plane sums perpendicular to orient direction
        if ( orient == 1 ) p(i,:) = mean(mean(Y1,3),2); end
        if ( orient == 2 ) p(i,:) = mean(mean(Y1,1),3); end
        if ( orient == 3 ) p(i,:) = mean(mean(Y1,1),2); end
    end
    % Select a low value for each plane, and set threshold a bit above it.
    pq = 0.5*( min(p) + median(p,1));
    % Preview estimate of bad slice fraction...
    prebad = length(find(p(1,:) > pq + OUTSLICEdef));
    prebad = length(find(p(2,:) > pq + OUTSLICEdef)) + prebad;
    prebad = length(find(p(3,:) > pq + OUTSLICEdef)) + prebad;
    percentbad = round(prebad*100/(3*nslice));
    disp('Estimated percent bad slices at default threshold')
    disp(percentbad)
    % User Input Threshold, and default suggestion.
    if exist('OUTSLICE','var') == 0
        OUTSLICE = spm_input(['Select threshold (default is shown)' ],1,'n',OUTSLICEdef);
    end
    pq = pq + OUTSLICE;  % pq array is the threshold test for bad slices.
    fprintf('\n Interpolating new values for bad slices when ');
    fprintf('\n average value outside head is %3d counts over baseline.\n',OUTSLICE);
    fprintf('\n Bad slice log will be written.');
    fprintf('\n');
    %  DETECT, COUNT, AND LOG THE ARTIFACTS
    disp('Writing Artifact Log to location of image set')
    [ dirname, sname ] = fileparts(P(1).fname);   %  XXXXXX
    tstamp = clock;
    % % % % % %     filen = ['BadSliceLog',date,'Time',num2str(tstamp(4)),num2str(tstamp(5)),'.txt'];
    filen = ['artslice_log_',tmp_img1_name,'.txt'];
    logname = fullfile(dirname,filen);
    logID = fopen(logname,'wt');
    fprintf(logID,'Bad Slice Log, Image Set first file is:  \n  %s\n', sname);
    fprintf(logID,'Test Criteria for Artifacts');
    fprintf(logID,'\n  Outslice (counts) = %4.1f' , OUTSLICE);
    fprintf(logID,'\n Slice Artifact List (Vol 1 = First File)\n');
    BadError = 0;     % Counts bad slices
end

if allfilt == 1
    fprintf('\n MEDIAN FILTER ALL THE DATA');
    fprintf('\n All scans are being filtered by median 3-point time filter.');
    fprintf('\n This is safe to do when TR = 2. Otherwise, use the');
    fprintf('\n    Remove Bad Slices option instead.\n');
    fprintf('\n');
end
spm_input('!DeleteInputObj');


%% Main Loop - Filter everything but the first and last volumes.
% Procedure:
% first two volumes in Y4(1) and Y4(2)
% start 'for loop' from third volume and put it in Y4(3)
% mask Y4(2) with automask and put in Y
% find mean(Y) and compare with threshold for specific slice pq(j)
% if higher than threshold interpolate slice as mean from adjacent Y4(1) and Y4(3)
% shift Y4(1<-2<-3) and check next Y4(3)


nscans = size(P,1);
STDvol = shi_art_slice_STD(P, orient, nscans, pq, Automask, vx, vy, vz);
Y4(1,:,:,:) = spm_read_vols(P(1));  % rows vary fastest
Y4(2,:,:,:) = spm_read_vols(P(2));
for i = 3:nscans
    Y4(3,:,:,:) = spm_read_vols(P(i));  % this is the post volume
    if allfilt == 1  % Filter all data by median.
        Yn = median(Y4,1);
        Yn2 = squeeze(Yn);
    end
    if slicefilt == 1  % Repair bad slices by linear interpolation.
        % Check if outside head value is over pq. If so, filter.

        Yn1 = squeeze(Y4(1,:,:,:)); % pre-post volumes
        Yn3 = squeeze(Y4(3,:,:,:));
        Ypre = ( 1 - Automask).*Yn1;
        Ypost = ( 1 - Automask).*Yn3;

        Yn2 = squeeze(Y4(2,:,:,:)); % actual volume
        Y = ( 1 - Automask).*Yn2;
        if ( orient == 1 )
            pypre = mean(mean(Ypre,3),2);
            pypost = mean(mean(Ypost,3),2);
            py = mean(mean(Y,3),2);
            for j = 1:vx
                if py(j) > pq(j)  % slice needs to be filtered.
                    % new code, check pre post slices are not over threshold
                    if pypost(j) > pq(j) && pypre(j) < pq(j)      % POST BAD, PRE OK

                        % calculate amount of change this slice would
                        % have changed according to the change of the
                        % adjacent slice
                        addvalue = 0;
                        if j > 1                                                    % if slice not first
                            change = LastGoodVol(j-1,:,:) - Yn2(j-1,:,:);           % how much the adjacent slice changed between last volume and this volume
                            sliceSTD = STDvol(j-1,:,:);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(j,:,:).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        % if slice is first, addvalue =0  and
                        % slice would be equal to pre in time
                        Yn2(j,:,:) = squeeze(Y4(1,j,:,:))+addvalue;       % set to pre + calculated change from adjacent slice
                        fprintf('   Interpolated sagittal with pre. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated sagittal with pre. Vol %d, Slice %d.\n',i-1,j);

                    elseif pypost(j) < pq(j) && pypre(j) > pq(j)   % POST OK, PRE BAD

                        % calculate amount of change this slice would
                        % have changed according to the adjacent slice
                        addvalue = 0;
                        if j > 1 && pypost(j-1) < pq(j-1)                           % if slice not first and adjacent of next volume not bad
                            change = squeeze(Y4(3,j-1,:,:)) - Yn2(j-1,:,:);         % how much the adjacent slice will change to next volume
                            sliceSTD = STDvol(j-1,:,:);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(j,:,:).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        % if slice is first or adjacent of next volume is bad
                        % addvalue=0 and slice will be equal to post in time
                        Yn2(j,:,:) = squeeze(Y4(3,j,:,:)) + addvalue;        % set to post + calculated change from adjacent slice
                        fprintf('   Interpolated sagittal with post. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated sagittal with post. Vol %d, Slice %d.\n',i-1,j);

                    elseif pypost(j) > pq(j) && pypre(j) > pq(j)   % PRE BAD, POST BAD

                        % calculate amount of change this slice would
                        % have changed according to the change of the
                        % adjacent slice
                        addvalue = 0;
                        if j > 1                                                    % if slice not first
                            change = LastGoodVol(j-1,:,:) - Yn2(j-1,:,:);           % how much the adjacent slice changed between last volume and this volume
                            sliceSTD = STDvol(j-1,:,:);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(j-1,:,:).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        Yn2(j,:,:) = LastGoodVol(j,:,:) + addvalue;          % set to last good known slice + value of change
                        fprintf('   Interpolated sagittal with latest good one. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated sagittal with latest good one. Vol %d, Slice %d.\n',i-1,j);

                    else                                                        % PRE GOOD, POST GOOD (original pre-post inrepolation)

                        Yn2(j,:,:) = squeeze((Y4(1,j,:,:) + Y4(3,j,:,:))/2.0);  % interpolate pre and post
                        fprintf('   Interpolated sagittal. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated sagittal. Vol %d, Slice %d.\n',i-1,j);

                    end

                    LastGoodVol = Yn2(j,:,:);
                    BadError = BadError + 1;
                end
                LastGoodVol = Yn2;
            end
        end
        if ( orient == 2 )
            pypre = mean(mean(Ypre,1),3);
            pypost = mean(mean(Ypost,1),3);
            py = mean(mean(Y,1),3);
            for j = 1:vy
                if py(j) > pq(j)  % slice needs to be filtered.
                    % new code, check pre post slices are not over threshold
                    if pypost(j) > pq(j) && pypre(j) < pq(j)      % post bad, pre ok

                        % calculate amount of change this slice would
                        % have changed according to the change of the
                        % adjacent slice
                        addvalue = 0;
                        if j > 1                                                    % if slice not first
                            change = LastGoodVol(:,j-1,:) - Yn2(:,j-1,:);           % how much the adjacent slice changed between last volume and this volume
                            sliceSTD = STDvol(:,j-1,:);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(:,j,:).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        % if slice is first, addvalue =0  and
                        % slice would be equal to pre in time
                        Yn2(:,j,:) = squeeze(Y4(1,:,j,:))+addvalue;       % set to pre + calculated change from adjacent slice
                        fprintf('   Interpolated coronal with pre. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated coronal with pre. Vol %d, Slice %d.\n',i-1,j);
                    elseif pypost(j) < pq(j) && pypre(j) > pq(j)   % post ok, pre bad

                        % calculate amount of change this slice would
                        % have changed according to the adjacent slice
                        addvalue = 0;
                        if j > 1 && pypost(j-1) < pq(j-1)                           % if slice not first and adjacent of next volume not bad
                            change = squeeze(Y4(3,:,j-1,:)) - Yn2(:,j-1,:);         % how much the adjacent slice will change to next volume
                            sliceSTD = STDvol(:,j-1,:);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(:,j,:).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        % if slice is first or adjacent of next volume is bad
                        % addvalue=0 and slice will be equal to post in time
                        Yn2(:,j,:) = squeeze(Y4(3,:,j,:)) + addvalue;        % set to post + calculated change from adjacent slice
                        fprintf('   Interpolated coronal with post. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated coronal with post. Vol %d, Slice %d.\n',i-1,j);
                    elseif pypost(j) > pq(j) && pypre(j) > pq(j)   % pre bad, post bad

                        % calculate amount of change this slice would
                        % have changed according to the change of the
                        % adjacent slice
                        addvalue = 0;
                        if j > 1                                                    % if slice not first
                            change = LastGoodVol(:,j-1,:) - Yn2(:,j-1,:);           % how much the adjacent slice changed between last volume and this volume
                            sliceSTD = STDvol(:,j-1,:);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(:,j-1,:).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        Yn2(:,j,:) = LastGoodVol(:,j,:) + addvalue;          % set to last good known slice + value of change
                        fprintf('   Interpolated coronal with latest good one. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated coronal with latest good one. Vol %d, Slice %d.\n',i-1,j);
                    else                                                        % pre good, post good
                        Yn2(:,j,:) = squeeze((Y4(1,:,j,:) + Y4(3,:,j,:))/2.0);  % interpolate pre and post
                        fprintf('   Interpolated coronal. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated coronal. Vol %d, Slice %d.\n',i-1,j);
                    end

                    LastGoodVol = Yn2(:,j,:);
                    BadError = BadError + 1;
                end
            end
        end
        if ( orient == 3 )
            pypre = mean(mean(Ypre,1),2);
            pypost = mean(mean(Ypost,1),2);
            py = mean(mean(Y,1),2);
            for j = 1:vz
                if py(j) > pq(j)  % slice needs to be filtered.
                    % new code, check pre post slices are not over
                    % threshold and interpolates intelligently

                    if pypost(j) > pq(j) && pypre(j) < pq(j)                         % post bad, pre ok

                        % calculate amount of change this slice would
                        % have changed according to the change of the
                        % adjacent slice
                        addvalue = 0;
                        if j > 1                                                    % if slice not first
                            change = LastGoodVol(:,:,j-1) - Yn2(:,:,j-1);           % how much the adjacent slice changed between last volume and this volume
                            sliceSTD = STDvol(:,:,j-1);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(:,:,j).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        % if slice is first, addvalue =0  and
                        % slice would be equal to pre in time
                        Yn2(:,:,j) = squeeze(Y4(1,:,:,j)) + addvalue;                 % set to pre + calculated change from adjacent slice
                        fprintf('   Interpolated axial with pre. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated axial with pre. Vol %d, Slice %d.\n',i-1,j);


                    elseif pypost(j) < pq(j) && pypre(j) > pq(j)   % post ok, pre bad

                        % calculate amount of change this slice would
                        % have changed according to the adjacent slice
                        addvalue = 0;
                        if j > 1 && pypost(j-1) < pq(j-1)                           % if slice not first and adjacent of next volume not bad
                            change = squeeze(Y4(3,:,:,j-1)) - Yn2(:,:,j-1);         % how much the adjacent slice will change to next volume
                            sliceSTD = STDvol(:,:,j-1);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(:,:,j).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        % if slice is first or adjacent of next volume is bad
                        % addvalue=0 and slice will be equal to post in time
                        Yn2(:,:,j) = squeeze(Y4(3,:,:,j)) + addvalue;        % set to post + calculated change from adjacent slice
                        fprintf('   Interpolated axial with post. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated axial with post. Vol %d, Slice %d.\n',i-1,j);


                    elseif pypost(j) > pq(j) && pypre(j) > pq(j)   % pre bad, post bad

                        % calculate amount of change this slice would
                        % have changed according to the change of the
                        % adjacent slice
                        addvalue = 0;
                        if j > 1                                                    % if slice not first
                            change = LastGoodVol(:,:,j-1) - Yn2(:,:,j-1);           % how much the adjacent slice changed between last volume and this volume
                            sliceSTD = STDvol(:,:,j-1);                             % what is the STD for the adjacent slice
                            STDchange = change./sliceSTD;                           % ratio of change of adjacent slice compared to its STD
                            addvalue = STDvol(:,:,j).*STDchange;                    % apply the same ratio to calculate current slice change
                        end

                        Yn2(:,:,j) = LastGoodVol(:,:,j) + addvalue;          % set to last good known slice + value of change
                        fprintf('   Interpolated axial with latest good one. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Replaced axial with latest good one. Vol %d, Slice %d.\n',i-1,j);

                    else                                                        % pre good, post good
                        Yn2(:,:,j) = squeeze((Y4(1,:,:,j) + Y4(3,:,:,j))/2.0);  % interpolate pre and post
                        fprintf('   Interpolated axial. Vol %d, Slice %d.\n',i-1,j);
                        fprintf(logID,'   Interpolated axial. Vol %d, Slice %d.\n',i-1,j);

                    end
                    BadError = BadError + 1;
                end
            end
        end
    end

    LastGoodVol = Yn2;

    if  ( allfilt == 0 & slicefilt == 0 )  % Head mask only.
        Yn2 = squeeze(Y4(2,:,:,:));
    end
    % Yn is a 4D file including the smoothed Y2 now.
    if ( amask ) Yn2 = Yn2.*Automask; end
    % Prepare the header for the filtered volume.
    V = spm_vol(P(i-1).fname);
    v = V;
    [dirname, sname, sext ] = fileparts(V.fname);
    sfname = [ prechar, sname ];
    filtname = fullfile(dirname,[sfname sext]);
    v.fname = filtname;
    %spm_write_vol(v,Yn2);
    noscale_write_vol(v,Yn2);
    showprog = [' Writing volume   ', prechar, sname, sext ];
% % % % % % % %     disp(showprog);
    if i == 3   % Write unfiltered first scan
        Yn1 = squeeze(Y4(1,:,:,:));
        if ( amask ) Yn1 = Yn1.*Automask; end
        [dirname, sname, sext ] = fileparts(P(1).fname);
        sfname = [ prechar, sname ];
        filtname = fullfile(dirname,[sfname sext]);
        v.fname = filtname;
        %spm_write_vol(v,Yn1);
        noscale_write_vol(v,Yn1);
        showprog = [' Writing volume   ', prechar, sname, sext ];
        disp(showprog);
    end
    if i == nscans  % Write unfiltered last scan.
        Yn1 = squeeze(Y4(3,:,:,:));
        if ( amask ) Yn1 = Yn1.*Automask; end
        Vlast = spm_vol(P(nscans).fname);
        [dirname, sname, sext ] = fileparts(P(nscans).fname);
        sfname = [ prechar, sname ];
        filtname = fullfile(dirname,[sfname sext]);
        v.fname = filtname;
        %spm_write_vol(v,Yn1);
        noscale_write_vol(v,Yn1);
        showprog = [' Writing volume   ', prechar, sname, sext ];
        disp(showprog);
    end

    % Slide the read volumes window up.
    Y4(1,:,:,:) = Y4(2,:,:,:);
    Y4(2,:,:,:) = Y4(3,:,:,:);
end

%  Summarize the slice errors, in the bad slice case.
if slicefilt == 1
    totalslices = nscans*nslice;
    ArtifactCount = BadError;  %  Number of unique bad slices with artifacts
    badpercent = BadError*100.0/totalslices;
    if ( ArtifactCount == 0 )
        disp(' CLEAN DATA! No bad slices detected');
        fprintf(logID, '\n\nCLEAN DATA. No bad slices detected.');
    end
    if  ArtifactCount > 0
        fprintf(logID,'\n\n Number of slices repaired = %4.0d',ArtifactCount);
        fprintf(logID,'\n\n Percentage of slices repaired = %5.1f',badpercent);
    end
    fclose(logID);

    report_fname = fullfile(dirname,['artslice_summ_',tmp_img1_name,'.txt']);
    fid = fopen(report_fname,'wt');
    fprintf(fid,'slice_total = %g\n',totalslices);
    fprintf(fid,'slice_bad   = %g\n',ArtifactCount);
    fprintf(fid,'percent_bad = %g\n',badpercent);
    fclose(fid);
end

disp(['Done! ' num2str(BadError) ' slices fixed (' num2str(round(badpercent)) '%)' ]);



%----------------------------------

function  Y = shi_art_automask(Image,Threshold,WriteVol)
% Y = art_automask( Image, Threshold, WriteVol )      (v2.3)
% art_automask;
%
%    Calculates a pretty good mask image from an input full-head Image
% in order to detect artifacts in the data. The threshold
% is higher than usual for SPM because the spiral scan generates more
% noise outside the head. When called with one
% or no arguments, the mask is adapted for the input image, and the mask is
% written to "ArtifactMask.img" for later user review. The adaptation
% sets a threshold level, removes speckle, and removes corner artifacts.
% The adaptive mask is usually a bit larger than the head.
%    The GUI version allows a user to select a threshold,e.g. by looking at
% an SPM display of the image to estimate a value. The program will fill
% in small holes in the mask, so you can pick slightly high values to
% suppress noise and still obtain a mask without gaps.
%
% GUI INPUT
%  Image  - select the image to use
%  Threshold  - Select a threshold to use.
%  UserMask  - Select a name for the derived mask file.
% BATCH INPUT:
%  Image  - full path name of an image file, e.g. 'C:\test\V002.img'
%  Threshold- measured as a FRACTION of the range of Image
%      If > 0, applies that threshold. Values from 0.10 to 0.30 are typical.
%         e.g. if the Image range is [0,2000] and Threshold = 0.15,
%         then a fixed threshold of 300 is applied.
%      If not supplied, threshold is set adaptively to a value that is
%         0.2 to 0.4 of the max range of the smoothed image.
%  WriteVol = 0 does not write a mask file.
%      If not supplied, mask file is written.
%OUTPUT:
%  Y  3D mask array, with 1's and 0's.
%  ArtifactMask.img, written to the image directory, if WriteVol = 1.
%  UserMaskname.img, written to the image directory in GUI maode.
%
% Paul Mazaika  April 2004.
% V.2  adapts to different mean image levels.  Paul Mazaika Aug. 2006.
% V2.1 write maskvalue=1. get/write changes for SPM5. (2/07)
% v2.3 support SPM12 (12/14)


% Adaptive Threshold Logic depends on estimated error rate outside the head.
% Fraction of points outside the head that pass the threshold must be small.
% For each slice, set the slice to zero if the fraction of mask points on
% the slice is smaller than parameter FMR.
FMR = 0.015;   % False Mask Percentage.

% Configure while preserving old SPM versions
spmv = spm('Ver'); spm_ver = 'spm5';  % chooses spm_select to read vols
if (strcmp(spmv,'SPM2')) spm_ver = 'spm2'; end
if (strcmp(spmv,'SPM2') || strcmp(spmv,'SPM5')) spm_defaults;
else spm('Defaults','fmri'); end


% Get the image data.
V = spm_vol(Image(:,:));  % Input could be a matrix; only need one image.
n  = prod(size(V));
Y = spm_read_vols(V);
% Fill in the small holes and reduce the noise spikes.
Y = smooth3(Y);  % default 3x3x3 box smoothing.
Yr = max(max(max(Y))) - min(min(min(Y)));  % previously used range.

% User defined mask threshold
if Threshold > 0  % Make the mask directly
    % Array temp is logical array with 1's and 0's
    temp(:,:,:) = (Y(:,:,:)>Threshold);
end

% Adaptive Mask Threshold
if ( Threshold == -1 )   % Find a threshold that discards three outer faces.
    % Use gray matter density as lower limit- count is 400.
    Tlow = fix(0.2*Yr);  Thigh = fix(0.4*Yr); Tskip = max(fix(Tlow/20),1); % upper thresh 0.5 was 0.4 before
    for Tbar = Tlow:Tskip:Thigh   % 400:20:800
        temp(:,:,:) = (Y(:,:,:) > Tbar);
    	% Count the number of mask points in the far faces of the volume
        xdim = size(Y);
        count1 = sum(sum(temp(:,:,1)));
        count2 = sum(sum(temp(:,:,xdim(3))));
        count3 = sum(sum(temp(:,1,:)));
        count4 = sum(sum(temp(:,xdim(2),:)));
        count5 = sum(sum(temp(1,:,:)));
        count6 = sum(sum(temp(xdim(1),:,:)));
        % Always have one face with large counts, sometimes have 2 such faces.
        countA = count1+count2+count3+count4+count5+count6;
        Xbar = [ count1 count2 count3 count4 count5 count6 ];
        Ybar = sort(Xbar);
        countC = Ybar(1) + Ybar(2);  % the two smallest face counts
        countB = Ybar(1) + Ybar(2) + Ybar(3);  % three smallest face counts
        % Number of voxels on 3 faces is approximately:
        nvox = xdim(1)*xdim(2) + xdim(2)*xdim(3) + xdim(1)*xdim(3);
        if ( countC < FMR*nvox )
            break;   % Exit the For loop, current Tbar is good enough.
        end
        %          iter = iter + 1;
        %          ygraph(1,iter) = countA;
        %          ygraph(2,iter) = countB;
        %          ygraph(3,iter) = countC;
    end
    disp('Adaptive Mask Threshold')
    disp(Tbar)
    %disp(countC)
    %disp(ygraph')
end
if Threshold == -1
    if Tbar >= Thigh-Tskip
        disp('Automask program failed.   Try choosing a mean image,')
        disp(' or manually set a threshold. Type help art_automask.')
        return
    end
end

% Clean up the corner alias artifact sometimes evident in the spiral data.
% If the edge of a plane has much more data than the center, delete the
% edge. Then if any plane has too few data points (meaning close to noise level)
% then set the entire plane in the mask image to zero. Note for spiral scans
% a better idea might be to determine the orientation of the spiral, and just
% remove the potential alias effects at those four edges corners.
xdim = size(Y);
iedge = floor(max(5,0.2*xdim(1)));
jedge = floor(max(5,0.2*xdim(2)));
kedge = floor(max(5,0.2*xdim(3)));
% Clear out the edges, and then check the planes.
for i = 1:xdim(1)
    % If edges are bigger than center, then delete edges.
    fmaski = sum(sum(temp(i,jedge:xdim(2)-jedge,kedge:xdim(3)-kedge)));
    mo1 = sum(sum(temp(i,1:jedge,:)));
    if ( mo1 > 2*fmaski ) temp(i,1:5,:) = 0; end
    mo2 = sum(sum(temp(i,xdim(2)-jedge:xdim(2),:)));
    if (mo2 > 2*fmaski) temp(i,xdim(2)-4:xdim(2),:) = 0; end
    mo3 = sum(sum(temp(i,:,1:kedge)));
    if (mo3 > 2*fmaski)  temp(i,:,1:5) = 0; end
    mo4 = sum(sum(temp(i,:,xdim(3)-kedge:xdim(3))));
    if (mo4 > 2*fmaski) temp(i,:,xdim(3)-4:xdim(3)) = 0; end
    % If face is about the noise level, then delete the face.
    fmaski = sum(sum(temp(i,:,:)));
    if fmaski < 2*FMR*xdim(2)*xdim(3)
        temp(i,:,:) = 0;
    end
end
for j = 1:xdim(2)
    fmaskj = sum(sum(temp(iedge:xdim(1)-iedge,j,kedge:xdim(3)-kedge)));
    mo1 = sum(sum(temp(1:iedge,j,:)));
    if ( mo1 > 2*fmaskj ) temp(1:5,j,:) = 0; end
    mo2 = sum(sum(temp(xdim(1)-iedge:xdim(1),j,:)));
    if (mo2 > 2*fmaskj) temp(xdim(1)-4:xdim(1),j,:) = 0; end
    mo3 = sum(sum(temp(:,j,1:kedge)));
    if (mo3> 2*fmaskj)  temp(:,j,1:5) = 0; end
    mo4 = sum(sum(temp(:,j,xdim(3)-kedge:xdim(3))));
    if (mo4 > 2*fmaskj) temp(:,j,xdim(3)-4:xdim(3)) = 0; end
    fmaskj = sum(sum(temp(:,j,:)));
    if fmaskj < 2*FMR*xdim(1)*xdim(3)
        temp(i,:,:) = 0;
    end
end
for k = 1:xdim(3)
    fmaskk = sum(sum(temp(:,:,k)));
    if fmaskk < 2*FMR*xdim(2)*xdim(1)
        temp(i,:,:) = 0;
    end
end

% Outputs
Y = temp;

if ( WriteVol == 1 )
    v = V;  % preserves the header structure
    [dirname, xname, xext ] = fileparts(V.fname);
    artifname = ['artslice_mask_', xname,  xext];
    artifpath = fullfile(dirname,artifname);
    v.fname = artifpath;
    noscale_write_vol(v,Y);
end




%---------------------------------------------------------------
% Create and write image without the scale and offset steps
% This function is spm_write_vol without error checking and scaling.
function noscale_write_vol(V,Y);
V = spm_create_vol(V);
for p=1:V.dim(3),
    V = spm_write_plane(V,Y(:,:,p),p);
end;
%V = spm_close_vol(V);  % not for SPM5








function STDvol = shi_art_slice_STD(P, orient, nscans, pq, Automask, vx, vy, vz)
% This function takes a list of filenames P and calculates a volume of the
% standard deviation of the timeseries without taking in consideration bad
% slices during calculation.
% D. Postina - 2012

for i = 1:nscans
    thisvol = spm_read_vols(P(i));
    Yn2 = squeeze(thisvol);
    Y = ( 1 - Automask).*Yn2;
    
    if ( orient == 1 )
        py = mean(mean(Y,3),2);
        for j = 1:vx
            if py(j) > pq(j)    % drop slice from variance calculation
                allvols(i,j,1:size(Yn2,2),1:size(Yn2,3)) = NaN;
            else
                allvols(i,j,1:size(Yn2,2),1:size(Yn2,3)) = Yn2(j,:,:);
            end
        end
    end;
    
    if ( orient == 2 )
        py = mean(mean(Y,1),3);
        for j = 1:vy
            if py(j) > pq(j)    % drop slice from variance calculation
                allvols(i,1:size(Yn2,1),j,1:size(Yn2,3)) = NaN;
            else
                allvols(i,1:size(Yn2,1),j,1:size(Yn2,3)) = Yn2(:,j,:);
            end
        end
    end;
    
    if ( orient == 3 )
        py = mean(mean(Y,1),2);
        for j = 1:vz
            if py(j) > pq(j)    % drop slice from variance calculation
                allvols(i,1:size(Yn2,1),1:size(Yn2,2),j) = NaN;
            else
                allvols(i,1:size(Yn2,1),1:size(Yn2,2),j) = Yn2(:,:,j);
            end
        end
    end;
end
temp = art_nanstd(allvols);
STDvol(:,:,:) = temp(1,:,:,:);
clear allvols;


function [fff_std] = art_nanstd(data);
%
%   [f_std] = nanstd(data);
%
%Function which calculates the std (not NaN) of data containing
%NaN's.  NaN's are excluded completely from calculation.

[m,n] = size(data);

for index = 1:n;
    not_nans = find(isnan(data(:,index)) == 0);
        if length(not_nans) > 0;
            f_std(index) = std(data(not_nans,index));
        else
            f_std(index) = NaN;
        end
end

% we have the 1xN array of STDs, lets put them in the volume
fff_std = zeros(1,size(data,2),size(data,3),size(data,4));
fff_std(:) = f_std(:);
