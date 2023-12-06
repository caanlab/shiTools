function imgout = shiSpmEffSize_CohensD_PairedSampleT(img1,img2,outname)

% calculates Cohen's d for difference between paired imaging datasets
%
% Zhenhao Shi 2020-4-5
%

n = numel(img1);
n2 = numel(img2);

if n ~= n2
    error('unmatched case number');
end

shiSpmImgCalc0([img1(:);img2(:)],sprintf('mean(X(1:%d,:)-X(%d:%d,:))./std(X(1:%d,:)-X(%d:%d,:))',n,n+1,2*n,n,n+1,2*n),outname);
imgout = outname;
