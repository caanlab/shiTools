function Motion = shiSpmPreprocMotCalc(Param6)

% calculates parameters relevant to head motion
%
% Motion = shiSpmPreprocMotCalc(Param6)
%
% input:
%            Param6: 6 columns of motion parameters from SPM preprocessing
%                    (3 translation params in mm, 3 rotation params in rad)
% output fields: 
%            Motion.Abs                 6 raw absolute motion param
%            Motion.Fd                  6 frame-wise displacement param
%            Motion.AbsSq               6 Abs param squared
%            Motion.FdSq                6 Fd param squared
%            Motion.Fd_Power            Fd by Power et al.
%            Motion.Fd_Trans_VanDijk    Fd (translation) by Van Dijk et al.
%            Motion.Fd_Rotat_VanDijk    Fd (rotation) by Van Dijk et al.
%            Motion.Fd_Jenkinson        Fd by Jenkinson et al. (FSL)
%            Motion.Abs_Power           Abs adapted from Power et al.
%            Motion.Abs_Trans_VanDijk   Abs (translation) adapted from Van Dijk et al.
%            Motion.Abs_Rotat_VanDijk   Abs (rotation) adapted from Van Dijk et al.
%            Motion.Abs_Jenkinson       Abs adapted from Jenkinson et al. (FSL)
%
% Zhenhao Shi, 2019/12/16
%

if ischar(Param6) && exist(Param6,'file')
    Param6 = importdata(Param6);
end

Abs = bsxfun(@minus,Param6,mean(Param6));
% translation (x, y, z), rotation (u, v, w), demeaned

Fd = [0,0,0,0,0,0;Param6(2:end,:) - Param6(1:end-1,:)];
% x(t) - x(t-1), ...

AbsSq = Param6.^2;
% x(t)^2, ...

FdSq = Fd.^2;
% (x(t) - x(t-1))^2, ... (i.e. dx, dy, dz, du, dv, dw)

Fd_Power = abs(Fd(:,1)) + abs(Fd(:,2)) + abs(Fd(:,3)) + abs(Fd(:,4)*50) + abs(Fd(:,5)*50) + abs(Fd(:,6)*50);
% |dx| + |dy| + |dz| + 50*|du| + 50*|dv| + 50*|dw|

Fd_Trans_VanDijk = sqrt(Fd(:,1).^2 + Fd(:,2).^2 + Fd(:,3).^2);
% sqrt( dx^2 + dy^2 + dz^2 )

Fd_Rotat_VanDijk = acos((cos(Fd(:,4)).*cos(Fd(:,5))+cos(Fd(:,4)).*cos(Fd(:,6))+cos(Fd(:,5)).*cos(Fd(:,6))+sin(Fd(:,4)).*sin(Fd(:,6)).*sin(Fd(:,5))-1)./2);
% arccos( ( cos(du)cos(dv) + cos(du)cos(dw) + cos(dv)cos(dw) + sin(du)sin(dv)sin(dw) -1 ) /2 )

Abs_Power = abs(Abs(:,1)) + abs(Abs(:,2)) + abs(Abs(:,3)) + abs(Abs(:,4)*50) + abs(Abs(:,5)*50) + abs(Abs(:,6)*50);
% |x| + |y| + |y| + 50*|u| + 50*|v| + 50*|w|

Abs_Trans_VanDijk = sqrt(Abs(:,1).^2 + Abs(:,2).^2 + Abs(:,3).^2);
% sqrt( x^2 + y^2 + y^2 )

Abs_Rotat_VanDijk = acos((cos(Abs(:,4)).*cos(Abs(:,5))+cos(Abs(:,4)).*cos(Abs(:,6))+cos(Abs(:,5)).*cos(Abs(:,6))+sin(Abs(:,4)).*sin(Abs(:,6)).*sin(Abs(:,5))-1)./2);
% arccos( ( cos(u)cos(v) + cos(u)cos(w) + cos(v)cos(w) + sin(u)sin(v)sin(w) -1 ) /2 )


Fd_Jenkinson = zeros(size(Param6,1),1);
Abs_Jenkinson = zeros(size(Param6,1),1);

T0 = spm_matrix(Param6(1,:));
iT0 = pinv(T0);
xc = [0 0 0]';

for i = 1:size(Param6,1)
    
    T1 = spm_matrix(Param6(max(1,i-1),:));
    iT1 = pinv(T1);
    T2 = spm_matrix(Param6(i,:));

    M_rel = T2*iT1 - eye(4);
    A_rel = M_rel(1:3,1:3);
    t_rel = M_rel(1:3,4);
    Fd_Jenkinson(i) = sqrt(1/5*50^2*trace(A_rel.'*A_rel)+(t_rel+A_rel*xc(:)).'*(t_rel+A_rel*xc(:)));

    M_abs = T2*iT0 - eye(4);
    A_abs = M_abs(1:3,1:3);
    t_abs = M_abs(1:3,4);
    Abs_Jenkinson(i) = sqrt(1/5*50^2*trace(A_abs.'*A_abs)+(t_abs+A_abs*xc(:)).'*(t_abs+A_abs*xc(:)));
    
end

Motion.Abs                   =   Abs                ;
Motion.Fd                    =   Fd                 ;
Motion.AbsSq                 =   AbsSq              ;
Motion.FdSq                  =   FdSq               ;
Motion.Fd_Power              =   Fd_Power           ;
Motion.Fd_Trans_VanDijk      =   Fd_Trans_VanDijk   ;
Motion.Fd_Rotat_VanDijk      =   Fd_Rotat_VanDijk   ;
Motion.Fd_Jenkinson          =   Fd_Jenkinson       ;
Motion.Abs_Power             =   Abs_Power          ;
Motion.Abs_Trans_VanDijk     =   Abs_Trans_VanDijk  ;
Motion.Abs_Rotat_VanDijk     =   Abs_Rotat_VanDijk  ;
Motion.Abs_Jenkinson         =   Abs_Jenkinson      ;
