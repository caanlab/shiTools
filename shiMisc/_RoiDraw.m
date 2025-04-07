name = 'SPL';
Rnum = 198;
Lnum = 199;
R = shiSpmRoiDraw(struct('atlas','Neuromorphometrics','value',Rnum),['shi_Neuromorphometrics_R',name,'.nii']);
L = shiSpmRoiDraw(struct('atlas','Neuromorphometrics','value',Lnum),['shi_Neuromorphometrics_L',name,'.nii']);
shiSpmImgCalc0([R;L],'i1|i2',['shi_Neuromorphometrics_Bi',name,'.nii'],[],0);
