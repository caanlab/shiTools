function SliceOrder = shiDicomSliceOrder(dicomName,SliceNumber)

% reads a dicom file and tries to determine slice order
%
% SliceOrder = shiDicomSliceOrder(dicomName,SliceNumber)
%   returns the order of slice acquisition by reading the dicom file with
%   file name dicomName, given the number of slices SliceNumber
% 
% ###########
% Zhenhao Shi 2018-03-23
% ###########


infoD = dicominfo(dicomName);
try
    str = infoD.Private_0029_1020;
catch
    error('Cannot find field Private_0029_1020. Try shiDicomSliceTime');
end
xstr = char(str');
n = strfind(xstr, 'sSliceArray.ucMode');
[~, r] = strtok(xstr(n:n+100), '=');
ucmode = strtok(strtok(r, '='));
switch(ucmode)
    case '0x1'
        sliceorder = 'Ascending';
    case '0x2'
        sliceorder = 'Descending';
    case '0x4'
        sliceorder = 'Interleaved';
    otherwise
        sliceorder = 'Order undetermined';
end

switch(sliceorder)
    case 'Ascending'
        SliceOrder = 1:1:SliceNumber;
    case 'Descending'
        SliceOrder = SliceNumber:-1:1;
    case 'Interleaved'
        % Interleaved order depends on whether slice number is odd or even!
        if mod(SliceNumber,2)
            SliceOrder = [1:2:SliceNumber 2:2:SliceNumber];
        else
            SliceOrder = [2:2:SliceNumber 1:2:SliceNumber];
        end
        warning('CAREFUL! Ensure your interleaved order is correct! Try shiDicomSliceTime') %#ok<*WNTAG>
    otherwise
        error('BAD ORDER! Check your slice order, and/or set it manually! Try shiDicomSliceTime')
end
