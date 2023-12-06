function imgout = shiSpmEffSize_CohensD_OneSampleT(img,outname)

% calculates Cohen's d for difference between imaging data and zero
%
% Zhenhao Shi 2020-4-5
%

shiSpmImgCalc0(img,'mean(X)./std(X)',outname);
imgout = outname;

