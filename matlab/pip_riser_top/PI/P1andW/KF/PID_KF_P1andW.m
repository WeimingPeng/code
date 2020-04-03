% PI control of P1 on the nonlinear model and using S-Function
clc
clear all
warning ('off')


controller='PI';
estimator='KF';
feedback='';
place = 'C:\Users\tereseva\Dropbox\Master\TereseCode\PID\KF\P1andW5perc\results\';
state = 'state';
input= 'input';
output2= 'output';

%% **** System Inputs ****
par = v4_new_4d_parameters(); %All parameters defined in one file
z0 = 0.10; % NB!!! CHOKE OPENING
wG_in0 = 0.36;
wL_in0 = 8.64;
global u0;
u0 = [z0;wG_in0;wL_in0];

%% ******** Initialize the state variables *************
disp('Initializing...')
tic
[x0,y0,par] =  v4_new_4d_initialize(u0,par);
toc
disp('Initialization completed.');
disp('#########################');
%% ******** linear model for estimation and control *************
[A,B,C,D]=v4_new_4d_linmod_knut(x0,y0,u0,par); %analytical
C = C(1:4,:);
D = D(1:4,:);
sys =(ss(A,B,C,D));
set(sys,'inputname',{'Z','wG_in','wL_in'},'outputname',{'P1','P2','Q','W'},'statename',{'m_G1','m_L1','m_G2','m_L2'});
Gu = tf(sys({'P1','P2','Q','W'},'Z'));
Gd = tf(sys({'P1','P2','Q','W'},{'wG_in','wL_in'}));
%% ****Scaling****
%y output( P1,P2,W)
%d disturbance(wG_in, wL_in)
%u input (Z)
De=(diag([1 1 1 1])); 
Dd=(diag([0.03*wG_in0 wL_in0*wL_in0]));
Du=min(z0,1-z0);

Gu=(De^-1)*Gu*Du;
Gd=(De^-1)*Gd*Dd;
set(Gu,'inputname','Z','outputname',{'P1','P2','Q','W'});
set(Gd,'inputname',{'wG_in','wL_in'},'outputname',{'P1','P2','Q','W'});

%%
Cs = (De^-1)*C;
Cs1 = Cs(1,:);
Cs2 = Cs(2,:);
Cs3 = Cs(3,:);
Cs4 = Cs(4,:);
Ds1 = (De^-1)*D(:,1)*Du;
Ds2 = (De^-1)*D(:,2:3)*Dd;
Ds = [Ds1 Ds2];
Bs1 = B(:,1)*Du;
Bs2 = B(:,2:3)*Dd;
Bs = [Bs1 Bs2];
syse=ss(A,Bs,Cs,Ds,'inputname',{'Z','wG_in','wL_in'},'outputname',{'P1','P2','Q','W'},'statename',{'m_G1','m_L1','m_G2','m_L2'});
sysc1=minreal(ss(A,Bs1,Cs1,Ds1(1,:),'inputname','Z','outputname','P1','statename',{'m_G1','m_L1','m_G2','m_L2'}));

%% Defineing Variables
t_f = 4.5; %hours
tfinal=t_f*3600; %seconds
t0=0;
Ts = 1;
t = t0:Ts:tfinal-Ts;
nf = length(t);

load wG_wL_nm

n_wG_in = n_wG_in(1:nf);
n_wL_in = n_wL_in(1:nf);
n_m = [n_m(:,1:nf);n_m(3,1:nf)];

%% Kalman Filter 1 (Measurement: P1,W)
sensors=[1 4];
known=1;
n_d = [n_wG_in; n_wL_in];
Qn = cov(n_d'); % Covariance of process noises
Rn = cov(n_m(sensors,:)'); % Covariance of measurement noises
Nn = zeros(2,length(sensors));

[kest,L,P,M,Q]=kalman(syse,Qn,Rn,Nn,sensors,known);
%% ############ PID Controller ##############
Kp = 1;
Ki = 0.0001;
tic
sim('PID_KF_P1andW_mdl')
runningTime=toc

%% For plotting
y(1,:)=YandYhat(:,1)';
y_hat(1,:)=YandYhat(:,2)';
y(2,:)=P2';
y(3,:)=W(:,1)';
y_hat(3,:)=W(:,2)';
u_in=Z';






