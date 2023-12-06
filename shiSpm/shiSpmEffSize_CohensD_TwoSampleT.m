function imgout = shiSpmEffSize_CohensD_TwoSampleT(img1,img2,outname)

% calculates Cohen's d for difference between independent imaging datasets
%
% Zhenhao Shi 2020-4-5
%

n1 = numel(img1);
n2 = numel(img2);

shiSpmImgCalc0([img1(:);img2(:)],sprintf('(mean(X(1:%d,:))-mean(X(%d:%d,:)))./sqrt(((%d).*std(X(1:%d,:)).^2+(%d).*std(X(%d:%d,:)).^2)./(%d))',n1,n1+1,n1+n2,n1-1,n1,n2-1,n1+1,n1+n2,n1+n2-2),outname);
imgout = outname;
