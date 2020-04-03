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

sys =(ss(A,Bs,Cs,Ds));
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Gp','m_Lp','m_Gr','m_Lr'});
sysc = minreal(sys);
sysm = minreal(sys('P_wh','Z1'));
Gu = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z1'));
Gd = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},{'d1','d2'}));
Gd11=(Gd('P_wh','d1')); Gd12=(Gd('P_wh','d2'));
%% H-infinity ptimal controller for plant
Gc = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z1'));
%Gu = tf(sys('P1','Z'));
NMEAS = 1;
NCON = 1;

% Find H-infinity optimal controller

W1 = [Wp1 zeros(1,10)
      0 Wp2 zeros(1,9)
      zeros(1,2) Wp1 zeros(1,8)
      zeros(1,3) Wp1 zeros(1,7)
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
%P = P(1:34,:);
P = minreal(P([1:33 35],:));
% set(P,'InputGroup', struct('U1', 1:3, 'U2', 4));
% set(P,'OutputGroup',struct('Y1', 1:9, 'Y2', 10:12));
[K,CL,GAM,INFO] = hinfsyn(P,NMEAS,NCON,'TOLGAM',1e-6,'METHOD','ric','DISPLAY','on');

%[K,CL,GAM,INFO] = mixsyn(G,Wp,Wu,Wt);
disp(INFO)


disp('Transfer function of H-infinity controller for P_wh')
zpk(K)
L1 = sysm*(K);
S1 = 1/(1+L1);
T1 = 1 - S1;
R1 = K*S1;
KSGd11 = R1*Gd11;
KSGd12 = R1*Gd12;
SG1 = S1*sysm;
disp('Poles of sensitivity transfer function of H-infinity controller for P_wh')
pole(S1)


omega = logspace(-6,2,100);

[mag,phase] = bode(sysm,omega);
mag_sys = reshape(mag,length(omega),1);
phase_sys = reshape(phase,length(omega),1);

[mag,phase] = bode(K,omega);
mag_K = reshape(mag,length(omega),1);
phase_K = reshape(phase,length(omega),1);

[mag,phase] = bode(L1,omega);
mag_L1 = reshape(mag,length(omega),1);
phase_L1 = reshape(phase,length(omega),1);


allmargin(L1)


for i=1:length(omega)
    S1f(i) = abs(freqresp(S1,omega(i)));
    T1f(i) = abs(freqresp(T1,omega(i)));
    R1f(i) = abs(freqresp(R1,omega(i)));
    SG1_f(i) = abs(freqresp(SG1,omega(i)));
    KSGd11_f(i) = abs(freqresp(KSGd11,omega(i)));
    KSGd12_f(i) = abs(freqresp(KSGd12,omega(i)));
    Gd11f(i) = abs(freqresp(Gd11,omega(i)));
    Gd12f(i) = abs(freqresp(Gd12,omega(i)));
    Gd11invf(i) = 1/abs(freqresp(Gd11,omega(i)));
    Gd12invf(i) = 1/abs(freqresp(Gd12,omega(i)));
    Wp_inv(i) = abs(freqresp(GAM/Wp2,omega(i)));
    Wt_inv(i) = abs(freqresp(GAM/Wt,omega(i)));
    %Wu_inv(i) = abs(freqresp(GAM/Wu,omega(i)));
    %Wt_inv(i) = GAM/Wt;
    Wu_inv(i) = GAM/Wu;
end
disp('Peak of S:')
disp(max(S1f))
disp('Peak T:')
disp(max(T1f))
disp('Peak of KS:')
disp(max(R1f))
disp('Peak of KSGd1:')
disp(max(KSGd11_f))
disp('Peak of KSGd2:')
disp(max(KSGd12_f))
disp('Peak of SG:')
disp(max(SG1_f))
%%
figure(1)
clf
set(gcf,'Color',[1,1,1])
set(gca,'FontSize',8)

subplot(2,2,1)
set(gca,'FontSize',8,'OuterPosition',[0  0.5  0.5  0.5])
loglog(omega,mag_sys,'k',omega,Gd11f,'--r',omega,Gd12f,'--b','Linewidth',2);
title('Mag. of transfer function form Z_1 to P_{wh}','FontSize',10,'FontWeight','bold','FontName','Times');
legend('|G_1|','G_{d1}','G_{d2}')
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')
xlim([1e-6 1e2])

subplot(2,2,2)
set(gca,'FontSize',8,'OuterPosition',[0.5  0.5  0.5  0.5])
loglog(omega,mag_K,'k','Linewidth',2);
title('Mag. of controller gain form P_{wh} to Z_1','FontSize',10,'FontWeight','bold','FontName','Times');
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')
xlim([1e-6 1e2])

subplot(2,2,3)
set(gca,'FontSize',8,'OuterPosition',[0  0  0.5  0.5])
semilogx(omega,phase_sys,'k','Linewidth',2);
title('Phase. of transfer function form Z_1 to P_{wh}','FontSize',10,'FontWeight','bold','FontName','Times');
xlim([1e-6 1e2])
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')

subplot(2,2,4)
set(gca,'FontSize',8,'OuterPosition',[0.5  0  0.5  0.5])
semilogx(omega,phase_K,'k','Linewidth',2);
title('Phase of controller form P_{wh} to Z_1','FontSize',10,'FontWeight','bold','FontName','Times');
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')
xlim([1e-6 1e2])

%%
figure(2)
clf
set(gcf,'Color',[1,1,1])
set(gca,'FontSize',8)

subplot(2,2,1)
set(gca,'FontSize',8,'OuterPosition',[0  0.5  0.5  0.5])
loglog(omega,mag_L1,'k',omega,ones(size(omega)),'--r','Linewidth',2);
title('Mag. of the loop transfer function with P_{wh}','FontSize',10,'FontWeight','bold','FontName','Times');
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')
xlim([1e-6 1e2])

subplot(2,2,2)
set(gca,'FontSize',8,'OuterPosition',[0.5  0.5  0.5  0.5])
loglog(omega,S1f,'r',omega,Wp_inv,'--r',omega,T1f,'k',omega,Wt_inv,'--k','Linewidth',2)
legend('|S_{1}|','|\gamma/Wp|','|T_{1}|','|\gamma/Wt|')
title('Sensitivity and Complementary Sensitivity','FontSize',10,'FontWeight','bold','FontName','Times');
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')
xlim([1e-6 1e2])

subplot(2,2,3)
set(gca,'FontSize',8,'OuterPosition',[0  0  0.5  0.5])
semilogx(omega,phase_L1,'k',omega,180*ones(size(omega)),'--r','Linewidth',2);
title('Phase of the loop transfer function with P_{wh}','FontSize',10,'FontWeight','bold','FontName','Times');
xlim([1e-6 1e2])
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')

subplot(2,2,4)
set(gca,'FontSize',8,'OuterPosition',[0.5  0  0.5  0.5])
loglog(omega,R1f,'k',omega,Wu_inv,'--k',omega,KSGd11_f,'r',omega,KSGd12_f,'b','Linewidth',2);
title('Input usage KS','FontSize',10,'FontWeight','bold','FontName','Times');
legend('|KS_1|','|\gamma/Wu|','|KSGd_{11}|','|KSGd_{12}|')
xlabel('\omega [Rad/s]','FontSize',10,'FontName','Times')
xlim([1e-6 1e2])

%%
K = minreal(K);

d2_final = 1; % Disturbance 2 
sim('Hinf_Nonlinear_Z1_Pwh_sim')
%%
t_f = t(end);
figure(7)
clf
% rect = [0, 0, 14, 10];
% set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[14 10],'PaperPosition',rect)
% set(gca,'FontSize',8)
set(gcf,'Color',[1,1,1])
set(gca,'FontSize',16)
linewidth = 2;
labelsize = 16;
legsize = 20;

subplot(3,2,1)
set(gca,'FontSize',legsize,'OuterPosition',[0  2/3  0.5  1/3])
plot(t/3600,z1,'k',t/3600,z2,'r','Linewidth',linewidth);
legend('Z_1','Z_2',1,'Orientation','Vertical')
%title('a) Manipulated variable (choke valve opening Z)','FontSize',10,'FontWeight','bold','FontName','Times');
%ylabel('\itZ\rm [%]','FontSize',10,'FontName','Times')
%xlabel('time [h]','FontSize',10,'FontName','Times')
xlim([0 t_f/3600])
ylim([-0.5 0.5])

subplot(3,2,2)
set(gca,'FontSize',legsize,'OuterPosition',[0.5  2/3  0.5  1/3])
plot(t/3600,d1,'r',t/3600,d2,'k','Linewidth',linewidth);
legend('d_1', 'd_2',1,'Orientation','Vertical')
%title('Disturbances')
xlim([0 t_f/3600])
ylim([-1.1 1.1])

subplot(3,2,3)
set(gca,'FontSize',legsize,'OuterPosition',[0  1/3  0.5  1/3])
plot(t/3600,Win,'k','Linewidth',linewidth);
legend('W_{in}',1,'Orientation','Vertical')
xlim([0 t_f/3600])
ylim([-1 1])

subplot(3,2,4)
set(gca,'FontSize',legsize,'OuterPosition',[0.5  1/3  0.5  1/3])
plot(t/3600,Pbh,'k','Linewidth',linewidth);
legend('P_{bh}',1,'Orientation','Vertical')
xlim([0 t_f/3600])
ylim([-1 1])

subplot(3,2,5)
set(gca,'FontSize',legsize,'OuterPosition',[0  0  0.5  1/3])
plot(t/3600,rPwh,'r',t/3600,Pwh,'k','Linewidth',linewidth);
legend('set-point','P_{wh}',1,'Orientation','Vertical')
xlim([0 t_f/3600])
xlabel('Time [h]','FontSize',labelsize)
ylim([-1 1])

subplot(3,2,6)
set(gca,'FontSize',legsize,'OuterPosition',[0.5  0  0.5  1/3])
plot(t/3600,Prb,'k','Linewidth',linewidth);
legend('P_{rb}',1,'Orientation','Horizontal')
xlim([0 t_f/3600])
xlabel('Time [h]','FontSize',labelsize)
ylim([-1 1])