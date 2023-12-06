function SliceTime = shiDicomSliceTime(dicomName)

% reads a dicom file and tries to determine slice time (in millisecond)
%
% SliceTime = shiDicomSliceTime(dicomName)
%   returns the timing of slice acquisition by reading the dicom file with
%   file name dicomName. may only work for Siemens dicom data.
% 
% ###########
% Zhenhao Shi 2018-03-23
% ###########

hdr = spm_dicom_headers(dicomName);
SliceTime = hdr{1}.Private_0019_1029;