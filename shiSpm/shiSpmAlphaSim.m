function shiSpmAlphaSim(MaskFile, OutTxt, ConnectDef, SmoothFWHM, Voxel_P, n_Iteration)

% conducts conventional Monte Carlo simulation similar to the AlphaSim in AFNI (be aware of the Eklund paper!)
% 
% shiSpmAlphaSim(MaskFile, OutTxt, ConnectDef, SmoothFWHM, Voxel_P, n_Iteration)
%   Modified from http://www.restfmri.net
% 
%   Input:
%       MaskFile    - The image file indicate which voxels to analyze.
%       OutTxt      - The filename of the result file.
%       ConnectDef  - The number of neighborhood voxels that one voxel has.
%                      6, 18, or 18
%       SmoothFWHM  - Gaussian filter width (FWHM, in mm).
%       Voxel_P     - Individual voxel threshold probability
%       n_Iteration - Number of Monte Carlo simulations.
% 

VariableLine=10000;


VMask = spm_vol_nifti(MaskFile);
[Mask, XYZMask] = spm_read_vols(VMask);
VoxelDim = [
    (max(XYZMask(1,:)) - min(XYZMask(1,:)))/(length(unique(XYZMask(1,:)))-1);
    (max(XYZMask(2,:)) - min(XYZMask(2,:)))/(length(unique(XYZMask(2,:)))-1);
    (max(XYZMask(3,:)) - min(XYZMask(3,:)))/(length(unique(XYZMask(3,:)))-1);
    ]';
[nX,nY,nZ] = size(Mask);

Mask = logical(Mask);
nXYZ = numel(find(Mask));


if ConnectDef~=6 && ConnectDef~=18 && ConnectDef~=26
    error('connect must be 6, 18 or 26');
end
    
ft=zeros(1,VariableLine);   
mt=zeros(1,VariableLine);   
count=0;
suma=0;
sumsq=0;

fprintf('000000/%06d',n_Iteration);
for nt=1:n_Iteration
    foneimt=zeros(1,VariableLine); 
    fim=normrnd(0,1,nX,nY,nZ);
    fim = fim.*Mask;
    if SmoothFWHM ~= 0
      fim = gauss_filter(SmoothFWHM,fim,VoxelDim); 
    end
    fimca=reshape(fim,1,[]);
    count=count+nXYZ;
    suma=sum(fimca)+suma;
    sumsq=sum(fimca.*fimca)+sumsq;
    mean=suma/count;
    sd = sqrt((sumsq - (suma * suma)/count) / (count-1));
    
    zthr =norminv(1 - Voxel_P);
    xthr=sd*zthr+mean;
    fim(fim<=xthr)=0;      
    fim(fim>xthr)=1;
    a=numel(find(fim==1))/nXYZ;
    [theCluster, theCount] = bwlabeln(fim, ConnectDef);
    for i=1:theCount
        foneimt(numel(find(theCluster==i)))=foneimt(numel(find(theCluster==i)))+1;
    end
    for i=1:theCount
        ft(numel(find(theCluster==i)))=ft(numel(find(theCluster==i)))+1;
    end
    mt(find(foneimt,1,'last'))=mt(find(foneimt,1,'last'))+1;
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%06d/%06d',nt,n_Iteration);
end
g_max_cluster_size = find(mt, 1, 'last' );
total_num_clusters = sum(ft);
divisor=n_Iteration*nXYZ;
prob_table=zeros(1,g_max_cluster_size);
alpha_table=zeros(1,g_max_cluster_size);
cum_prop_table=zeros(1,g_max_cluster_size);
for i = 1:g_max_cluster_size
      prob_table(i) = i * ft(i) / divisor;
      alpha_table(i) = mt(i) / n_Iteration;
      cum_prop_table(i) = ft(i) / total_num_clusters;
end
for i = 1:g_max_cluster_size-1
      j = g_max_cluster_size - i +1;
      prob_table(j-1) = prob_table(j)+prob_table(j-1);
      alpha_table(j-1) = alpha_table(j)+alpha_table(j-1);
      cum_prop_table(i+1) = cum_prop_table(i)+cum_prop_table(i+1);
end


fid=fopen(sprintf('%s',OutTxt),'w');
 if(fid)
     if SmoothFWHM == 4.55
         SmoothFWHM=4;
     end
     
     fprintf(fid,'Mask filename = %s\n',MaskFile);
     fprintf(fid,'Voxels in mask = %d\n',nXYZ);
     fprintf(fid,'Gaussian filter width (FWHM, in mm) = %.3f\n',SmoothFWHM);
     fprintf(fid,'Cluster connection : connect = %.2f\n',ConnectDef);
     fprintf(fid,'Individual voxel threshold probability = %.3f\n',Voxel_P);
     fprintf(fid,'Number of Monte Carlo simulations = %d\n',n_Iteration);
     fprintf(fid,'Output filename = %s\n\n\n',OutTxt);
   
     fprintf(fid,'Cl Size\tFrequency\tCum Prop\tp/Voxel\tMax Freq\tAlpha\n');
     for i=1:g_max_cluster_size
         fprintf(fid,'%d\t\t%d\t\t%f\t%f\t%d\t\t%f\n',i,ft(i),cum_prop_table(i),prob_table(i),mt(i),alpha_table(i));
     end
     fclose(fid);
 end
 


function Q=gauss_filter(s,P,VOX)
if length(s) == 1; s = [s s s];end 
s  = s./VOX;					% voxel anisotropy
s  = max(s,ones(size(s)));			% lower bound on FWHM
s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

x  = round(6*s(1)); x = -x:x;
y  = round(6*s(2)); y = -y:y;
z  = round(6*s(3)); z = -z:z;
x  = exp(-(x).^2/(2*(s(1)).^2)); 
y  = exp(-(y).^2/(2*(s(2)).^2)); 
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
Q=P;
spm_conv_vol(P,Q,x,y,z,-[i,j,k]);
