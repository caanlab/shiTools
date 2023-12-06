function shiSpmDicomImport(DicomPatientDir, NiftiPatientDir, DicomFilePattern, NiftiSuffix)

% convers individual subject's dicom files to nifti files (cannot handle all subjects at once)
%
% shiSpmDicomImport(DicomPatientDir, NiftiPatientDir, DicomFilePattern, NiftiSuffix)
%   converts dicom files in path DicomPatientDir (string), filtered by
%   DicomFilePattern (string), to nifti files. Converted files are saved in
%   NiftiPatientDir (string). Dicom files will be searched for recursively
%   in DicomPatientDir and its subfolders. Default nifti suffix: .nii.
%
%   Example:
%       shiSpmDicomImport( ...
%           'D:\Study\Dicom\Sub01', ...
%           'D:\Study\Nifti\Sub01', ...
%           '*.dcm' );
%
%    ###########
% by Zhenhao Shi @ 2018-6-4
%    ###########
%

if ~exist('NiftiSuffix','var') || isempty(NiftiSuffix) || strcmpi(NiftiSuffix,'.nii') || strcmpi(NiftiSuffix,'nii')
    NiftiSuffix = 'nii';
elseif strcmpi(NiftiSuffix,'.img') || strcmpi(NiftiSuffix,'img')
    NiftiSuffix = 'img';
else
    error('NiftiSuffix should be either nii or img');
end

shiDisp({'Importing dicom files','','From:',['  ',DicomPatientDir],'To:',['  ',NiftiPatientDir]});
DicomPatientDir = char(shiFullFileName(DicomPatientDir));
NiftiPatientDir = shiMkdir(NiftiPatientDir);
cwd = pwd;

P = char(shiFileName_Recur(DicomPatientDir,DicomFilePattern));
hdr = spm_dicom_headers(P);

cd(NiftiPatientDir);
out = spm_dicom_convert(hdr,'all','series',NiftiSuffix);

outDir = cell(size(out.files));
for i = 1:numel(out.files)
    outDir{i} = fileparts(out.files{i});
end;
outDir = unique(outDir);

switch NiftiSuffix

    case 'nii'

        for i = 1:numel(outDir)
            fprintf('\n\n%s\n\n',outDir{i});
            cd(outDir{i});
            Nii = shiFileName('*.nii');
            NiiLeft = cell(size(Nii));
            for n = 1:numel(Nii)
                NiiLeft{n,1} = Nii{n}(1:18);
            end;
            assert(numel(unique(NiiLeft))==1,'unwanted nifti files found');
            assert(isequal(Nii,sort(Nii)),'nifti files badly sorted');
            for n = 1:numel(Nii)
                movefile(Nii{n},['vol_',num2str(n,'%.5d'),'.nii']);
                if mod(n,10)==1 || n==numel(Nii)
                    fprintf('%s <- %s\n',['vol_',num2str(n,'%.5d'),'.nii'],Nii{n});
                end
            end;
        end;

    case 'img'

        for i = 1:numel(outDir)
            fprintf('\n\n%s\n\n',outDir{i});
            cd(outDir{i});
            Img = shiFileName('*.img');
            Hdr = shiFileName('*.hdr');
            assert(numel(Img)==numel(Hdr),'.img and .hdr files unmatched');
            ImgLeft = cell(size(Img));
            HdrLeft = cell(size(Img));
            for n = 1:numel(Img)
                [~,imgName,~] = fileparts(Img{n});
                [~,hdrName,~] = fileparts(Hdr{n});
                assert(isequal(imgName,hdrName),'.img and .hdr files unmatched');
                ImgLeft{n,1} = Img{n}(1:18);
                HdrLeft{n,1} = Hdr{n}(1:18);
            end;
            assert(numel(unique(ImgLeft))==1,'unwanted nifti files found');
            assert(isequal(Img,sort(Img)),'nifti files badly sorted');
            for n = 1:numel(Img)
                movefile(Img{n},['vol_',num2str(n,'%.5d'),'.img']);
                movefile(Hdr{n},['vol_',num2str(n,'%.5d'),'.hdr']);
                if mod(n,10)==1 || n==numel(Img)
                    fprintf('%s <- %s\n',['vol_',num2str(n,'%.5d'),'.img'],Img{n});
                end
            end;
        end;

end;

cd(cwd);

