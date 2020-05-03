clc
clear all
warning ('off')
Current_Folder = cd('../WPR_model');
% **** System parameters ****
par = v2_new_5d_parameters(); 
% **** System Inputs ****
[z1,z2] = lin_point();
global u0;
u0 = [z1;z2;0;0];
% ******** Finding initial states *************
disp('Initializing...')
tic
[x0,y0,par] =  v2_new_5d_initialize(u0,par);
toc
disp('Initialization completed.');
disp('#########################');
%% ******** Linearization *************
[A,B,C,D]=v2_new_5d_linmodl_num(x0,u0,par);%numerical
sys =(ss(A,B,C,D));
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Gp','m_Lp','m_Gr','m_Lr'});
%% ****Scaling****
%y output( P1,P2,W)
%d disturbance(wG_in, wL_in)
%u input (Z)
[De,Dd,Du] = scaling(z1,z2);
[Wp1,Wp2,Wu,Wt] = Weights();
Bs1 = B(:,1:2)*Du;
Bs2 = B(:,3:4)*Dd;
Bs = [Bs1 Bs2];

Cs = (De^-1)*C;

Ds1 = (De^-1)*D(:,1:2)*Du;
Ds2 = (De^-1)*D(:,3:4)*Dd;
Ds = [Ds1 Ds2];

cd(Current_Folder);

sys =ss(A,Bs,Cs,Ds);
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Gp','m_Lp','m_Gr','m_Lr'});
sysc = minreal(sys);
Gu = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z2'));
Gd = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},{'d1','d2'}));
Gd231=Gd({'P_in','W_out'},'d1'); Gd232=Gd({'P_in','W_out'},'d2');
sysm = sysc({'P_in','W_out'},'Z2');
sysd = sysc({'P_in','W_out'},{'d1','d2'});
%% H-infinity ptimal controller for plant
Gc = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z2'));
%Gu = tf(sys('P1','Z'));
NMEAS = 2;
NCON = 1;

% Find H-infinity optimal controller

W1 = [Wp1 zeros(1,10)
      0 Wp1 zeros(1,9)
      zeros(1,2) Wp1 zeros(1,8)
      zeros(1,3) Wp2 zeros(1,7)
      zeros(1,4) Wp1 zeros(1,6)
      zeros(1,5) Wp1 zeros(1,5)
      zeros(1,6) Wp1 zeros(1,4)
      zeros(1,7) Wp1 zeros(1,3)
      zeros(1,8) Wp1 zeros(1,2)
      zeros(1,9) Wp1 0
      zeros(1,10) Wp1];
W2 = repmat(Wu,11,1);
W3 = [Wt zeros(1,10)
      0 Wt zeros(1,9)
      zeros(1,2) Wt zeros(1,8)
      zeros(1,3) Wt zeros(1,7)
      zeros(1,4) Wt zeros(1,6)
      zeros(1,5) Wt zeros(1,5)
      zeros(1,6) Wt zeros(1,4)
      zeros(1,7) Wt zeros(1,3)
      zeros(1,8) Wt zeros(1,2)
      zeros(1,9) Wt 0
      zeros(1,10) Wt];
%P = [W1 minreal(-W1*G); zeros(3,3) W2; zeros(3,3) W3*G ; eye(3) -G];
%P = minreal(ss(P));
P = augw(Gc,W1,W2,W3);
P = minreal(P([1:33 37 42],:));
% set(P,'InputGroup', struct('U1', 1:3, 'U2', 4));
% set(P,'OutputGroup',struct('Y1', 1:9, 'Y2', 10:12));
[K,CL,GAM,INFO] = hinfsyn(P,NMEAS,NCON,'TOLGAM',1e-6,'METHOD','ric','DISPLAY','on');
%[KK,CL,GAM,INFO] = mixsyn(G,Wp,Wu,Wt);
disp(GAM)
disp(INFO)

set(K,'inputname',{'P_in','W_out'},'outputname','Z2');
disp('Gamma:')
disp(GAM)
