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
sysm = minreal(sys({'P_wh','W_in'},'Z1'));
Gu = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z1'));
Gd = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},{'d1','d2'}));
Gd11=(Gd({'P_wh','W_in'},'d1')); Gd12=(Gd({'P_wh','W_in'},'d2'));
sysd = sysc({'P_wh','W_in'},{'d1','d2'});
%% H-infinity ptimal controller for plant
Gc = tf(sysc({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z1'));
%Gu = tf(sys('P1','Z'));
NMEAS = 2;
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
P = minreal(P([1:33 35 36],:));
% set(P,'InputGroup', struct('U1', 1:3, 'U2', 4));
% set(P,'OutputGroup',struct('Y1', 1:9, 'Y2', 10:12));
[K,CL,GAM,INFO] = hinfsyn(P,NMEAS,NCON,'TOLGAM',1e-6,'METHOD','ric','DISPLAY','on');

%[K,CL,GAM,INFO] = mixsyn(G,Wp,Wu,Wt);
disp(INFO)

set(K,'inputname',{'P_wh','W_in'},'outputname','Z1');
disp('Gamma:')
disp(GAM)

disp('Transfer function of H-infinity controller')
K_tf = tf(K);
zpk(K_tf)

omega = logspace(-6,1,100);
[mag,phase] = bode(K(1,1),omega);
mag2_dir = reshape(mag,length(omega),1);
phase2_dir = reshape(phase,length(omega),1);

[mag,phase] = bode(K(1,2),omega);
mag3_dir = reshape(mag,length(omega),1);
phase3_dir = reshape(phase,length(omega),1);

figure(1)
clf
subplot(2,2,1)
loglog(omega,mag2_dir,'k','Linewidth',2)
xlim([1e-6 1e1])
title('Mag form Pwh to Z')
%legend('With state estimation')
subplot(2,2,3)
semilogx(omega,phase2_dir,'k','Linewidth',2)
xlim([1e-6 1e1])
title('Phase form Pwh to Z')
%legend('With state estimation')
subplot(2,2,2)
loglog(omega,mag3_dir,'k','Linewidth',2)
xlim([1e-6 1e1])
title('Mag form W to Z')
%legend('With state estimation')
subplot(2,2,4)
semilogx(omega,phase3_dir,'k','Linewidth',2)
xlim([1e-6 1e1])
title('Phase form W to Z')
%legend('With state estimation')


G = sysm;
%************** Input sensitivity *************
LI=K*G;
SI = 1/(1+LI);
TI = 1-SI;
RI1 = K(1,1)*SI;
RI2 = K(1,2)*SI;
%*************** MIMO Sensitivity ****************
L=G*K;
S=eye(2)/(eye(2)+L);
T=eye(2)-S;
R = K*S;
KSGd = K*S*sysd;
SG = S*G;
%*********** Decentralized Sensitivity  *********
L11 = G(1,1)*K(1,1);
S11 = 1/(1+L11);
T11 = 1 - S11;
R11 = K(1,1)*S11;

L22 = G(2,1)*K(1,2);
S22 = 1/(1+L22);
T22 = 1 - S22;
R22 = K(1,2)*S22;

SG11 = S11*G(1,1);
SG22 = S22*G(2,1);

disp('Poles of S with H-infinity optimal controller for Pwh to Z')
pole(S(1,1))
disp('Poles of S with H-infinity optimal controller for W to Z')
pole(S(2,2))
%%
for i=1:length(omega)
    
    L11f(i) = abs(freqresp(L(1,1),omega(i)));
    S11f(i) = abs(freqresp(S(1,1),omega(i)));
    T11f(i) = abs(freqresp(T(1,1),omega(i)));
    R11f(i) = abs(freqresp(R(1,1),omega(i)));
    SG11f(i) = abs(freqresp(SG(1,1),omega(i)));
    KSGd1f(i) = abs(freqresp(KSGd(1,1),omega(i)));
    
    L12f(i) = abs(freqresp(L(1,2),omega(i)));
    S12f(i) = abs(freqresp(S(1,2),omega(i)));
    T12f(i) = abs(freqresp(T(1,2),omega(i)));
    R12f(i) = abs(freqresp(R(1,2),omega(i)));
    KSGd2f(i) = abs(freqresp(KSGd(1,2),omega(i)));
    
    L21f(i) = abs(freqresp(L(2,1),omega(i)));
    S21f(i) = abs(freqresp(S(2,1),omega(i)));
    T21f(i) = abs(freqresp(T(2,1),omega(i)));
    SG21f(i) = abs(freqresp(SG(2,1),omega(i)));
    
    L22f(i) = abs(freqresp(L(2,2),omega(i)));
    S22f(i) = abs(freqresp(S(2,2),omega(i)));
    T22f(i) = abs(freqresp(T(2,2),omega(i)));
    
    
    LIf(i) = abs(freqresp(LI,omega(i)));
    SIf(i) = abs(freqresp(SI,omega(i)));
    TIf(i) = abs(freqresp(TI,omega(i)));
    RI1f(i) = abs(freqresp(RI1,omega(i)));
    RI2f(i) = abs(freqresp(RI2,omega(i)));
    
    L11f_dec(i) = abs(freqresp(L11,omega(i)));
    S11f_dec(i) = abs(freqresp(S11,omega(i)));
    T11f_dec(i) = abs(freqresp(T11,omega(i)));
    R11f_dec(i) = abs(freqresp(R11,omega(i)));
    SG11f_dec(i) = abs(freqresp(SG11,omega(i)));
    
    L22f_dec(i) = abs(freqresp(L22,omega(i)));
    S22f_dec(i) = abs(freqresp(S22,omega(i)));
    T22f_dec(i) = abs(freqresp(T22,omega(i)));
    R22f_dec(i) = abs(freqresp(R22,omega(i)));
    SG22f_dec(i) = abs(freqresp(SG22,omega(i)));
    
    Wp_inv2(i) = abs(freqresp(GAM/Wp2,omega(i)));
    Wp_inv3(i) = GAM/Wp1; %abs(freqresp(GAM/Wp3,omega(i)));
    Wt_inv(i) = abs(freqresp(GAM/Wt,omega(i)));
    %Wu_inv(i) = abs(freqresp(GAM/Wu,omega(i)));
    %Wt_inv(i) = GAM/Wt;
    Wu_inv(i) = GAM/Wu;
    
end
disp('*********** Sensitivity of MIMO Transfer Functions ***********');
disp('Peak of MIMO S11:')
disp(max(S11f));

disp('Peak of MIMO T11:')
disp(max(T11f));

disp('Peak of MIMO S12:')
disp(max(S12f));

disp('Peak of MIMO T12:')
disp(max(T12f));

disp('Peak of MIMO S21:')
disp(max(S21f));

disp('Peak of MIMO T21:')
disp(max(T21f));

disp('Peak of MIMO S22:')
disp(max(S22f));

disp('Peak of MIMO T22:')
disp(max(T22f));

disp('Peak of MIMO KS for P_wh:')
disp(max(R11f));

disp('Peak of MIMO KS for W:')
disp(max(R12f));

disp('Peak of MIMO KSGd_1 (Fisrt Disturbance):')
disp(max(KSGd1f));

disp('Peak of MIMO KSGd_2 (Second Disturbance):')
disp(max(KSGd2f));

disp('Peak of MIMO SG11 for P_wh:')
disp(max(SG11f));

disp('Peak of MIMO SG21 for W:')
disp(max(SG21f));
disp('*********** Input Sensitivity ***********');
disp('Peak of MIMO SI (input sensitivity):')
disp(max(SIf));

disp('Peak of MIMO TI: (input compl-sensitivity)')
disp(max(TIf));

disp('Peak of MIMO KSI for P_wh:')
disp(max(RI1f));

disp('Peak of MIMO KSI for W:')
disp(max(RI2f));
disp('*********** Sensitivity of Decentralized Transfer Functions ***********');
disp('Peak of Decentralized S11:')
disp(max(S11f_dec));

disp('Peak of Decentralized T11:')
disp(max(T11f_dec));

disp('Peak of Decentralized S22:')
disp(max(S22f_dec));

disp('Peak of Decentralized T22:')
disp(max(T22f_dec));

disp('Peak of Decentralized SG11 for P_wh:')
disp(max(SG11f_dec));

disp('Peak of Decentralized SG22 for W:')
disp(max(SG22f_dec));
%%
figure(2)
clf
subplot(2,2,1)
loglog(omega,L11f,'k',omega,ones(size(omega)),'--r','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO Loop of P_wh')

subplot(2,2,2)
loglog(omega,L12f,'k',omega,ones(size(omega)),'--r','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO Off-diagonal Loop')

subplot(2,2,3)
loglog(omega,L21f,'k',omega,ones(size(omega)),'--r','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO off-digonal Loop')


subplot(2,2,4)
loglog(omega,L22f,'k',omega,ones(size(omega)),'--r','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO Loop of W')

%%
figure(3)
clf
subplot(2,2,1)
loglog(omega,S11f,'r',omega,Wp_inv2,'--r',omega,T11f,'k',omega,Wt_inv,'--k','Linewidth',2)
%loglog(omega,S11f,'r',omega,T11f,'k','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO Sensitivity and Complementary Sensitivity of P_wh')
legend('|S_{11}|','|\gamma/Wp|','|T_{11}|','|\gamma/Wt|')
%legend('|S_{11}|','|T_{11}|')

subplot(2,2,2)
loglog(omega,S12f,'r',omega,Wp_inv2,'--r',omega,T12f,'k',omega,Wt_inv,'--k','Linewidth',2)
%loglog(omega,S12f,'r',omega,T12f,'k','Linewidth',2)
xlim([1e-6 1e1])
title('off-digonal terms in MIMO Sensitivity and Complementary Sensitivity')
legend('|S_{12}|','|\gamma/Wp|','|T_{12}|','|\gamma/Wt|')
%legend('|S_{12}|','|T_{12}|')

subplot(2,2,3)
loglog(omega,S21f,'r',omega,Wp_inv3,'--r',omega,T21f,'k',omega,Wt_inv,'--k','Linewidth',2)
%loglog(omega,S21f,'r',omega,T21f,'k','Linewidth',2)
title('off-digonal terms in MIMO Sensitivity and Complementary Sensitivity')
xlim([1e-6 1e1])
legend('|S_{21}|','|\gamma/Wp|','|T_{21}|','|\gamma/Wt|')
%legend('|S_{21}|','|T_{21}|')

subplot(2,2,4)
loglog(omega,S22f,'r',omega,Wp_inv3,'--r',omega,T22f,'k',omega,Wt_inv,'--k','Linewidth',2)
%loglog(omega,S22f,'r',omega,T22f,'k','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO Sensitivity and Complementary Sensitivity of W')
legend('|S_{22}|','|\gamma/Wp|','|T_{22}|','|\gamma/Wt|')
%legend('|S_{22}|','|T_{22}|')

%%
figure(4)
clf
subplot(2,2,1)
loglog(omega,R11f,'k',omega,Wu_inv,'--k','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO KS of P_wh')
legend('|KS_{11}|','|\gamma/Wu|')

subplot(2,2,2)
loglog(omega,R12f,'k',omega,Wu_inv,'--k','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO KS of W')
legend('|KS_{12}|','|\gamma/Wu|')

subplot(2,2,3)
loglog(omega,KSGd1f,'k','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO KSGd_1 (Fisrt Disturbance)')
legend('|KSGd_1|')

subplot(2,2,4)
loglog(omega,KSGd2f,'k','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO KSGd_2 (Second Disturbance)')
legend('|KSGd_2|')
%%
figure(5)
clf
subplot(2,1,1)
loglog(omega,SIf,'r',omega,Wp_inv3,'--r',omega,TIf,'k',omega,Wt_inv,'--k','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO Input Sensitivity and Complementary Sensitivity')
legend('|\it\fontname{Times New Roman}S_I\rm|','|\gamma/Wp|','|\it\fontname{Times New Roman}T_I\rm|','|\gamma/Wt|')

subplot(2,1,2)
loglog(omega,RI1f,'k',omega,Wu_inv,'--k',omega,RI2f,'r','Linewidth',2)
xlim([1e-6 1e1])
title('MIMO Input KS(Input Usage)')
legend('|\it\fontname{Times New Roman}KS_{I1}\rm|','|\gamma/Wu|','|\it\fontname{Times New Roman}KS_{I2}\rm|')

%%
figure(6)
clf
subplot(2,2,1)
loglog(omega,S11f_dec,'r',omega,Wp_inv2,'--r',omega,T11f_dec,'k',omega,Wt_inv,'--k','Linewidth',2)
xlim([1e-6 1e1])
title('Decentralized Sensitivity and Complementary Sensitivity of P_wh')
legend('|S_{11}|','|\gamma/Wp|','|T_{11}|','|\gamma/Wt|')


subplot(2,2,2)
loglog(omega,S22f_dec,'r',omega,Wp_inv3,'--r',omega,T22f_dec,'k',omega,Wt_inv,'--k','Linewidth',2)
xlim([1e-6 1e1])
title('Decentralized Sensitivity and Complementary Sensitivity of W')
legend('|S_{22}|','|\gamma/Wp|','|T_{22}|','|\gamma/Wt|')

subplot(2,2,3)
loglog(omega,R11f_dec,'k',omega,Wu_inv,'--k','Linewidth',2)
xlim([1e-6 1e1])
title('Decentralized KS(Input Usage) of P_wh')
legend('|KS_{11}|','|\gamma/Wu|')

subplot(2,2,4)
loglog(omega,R22f_dec,'k',omega,Wu_inv,'--k','Linewidth',2)
xlim([1e-6 1e1])
title('Decentralized KS(Input Usage) of W')
legend('|KS_{22}|','|\gamma/Wu|')

%
K = minreal(K);
sim('Hinf_MISO_Linear_Z1_PwhWin_sim')
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
legsize = 12;

subplot(3,2,1)
set(gca,'FontSize',legsize,'OuterPosition',[0  2/3  0.5  1/3])
plot(t/3600,z1,'k',t/3600,z2,'r','Linewidth',linewidth);
legend('Z_1','Z_2',1,'Orientation','Vertical')
%title('a) Manipulated variable (choke valve opening Z)','FontSize',10,'FontWeight','bold','FontName','Times');
%ylabel('\itZ\rm [%]','FontSize',10,'FontName','Times')
%xlabel('time [h]','FontSize',10,'FontName','Times')
xlim([0 t_f/3600])
ylim([-1 1])

subplot(3,2,2)
set(gca,'FontSize',legsize,'OuterPosition',[0.5  2/3  0.5  1/3])
plot(t/3600,d1,'r',t/3600,d2,'k','Linewidth',linewidth);
legend('d_1', 'd_2',1,'Orientation','Vertical')
%title('Disturbances')
xlim([0 t_f/3600])
ylim([-1.1 1.1])

subplot(3,2,3)
set(gca,'FontSize',legsize,'OuterPosition',[0  1/3  0.5  1/3])
plot(t/3600,rWin,'r',t/3600,Win,'k','Linewidth',linewidth);
legend('set-point','W_{in}',1,'Orientation','Vertical')
xlim([0 t_f/3600])
ylim([-2 2])

subplot(3,2,4)
set(gca,'FontSize',legsize,'OuterPosition',[0.5  1/3  0.5  1/3])
plot(t/3600,Pbh,'k','Linewidth',linewidth);
legend('P_{bh}',1,'Orientation','Vertical')
xlim([0 t_f/3600])
ylim([-2 2])

subplot(3,2,5)
set(gca,'FontSize',legsize,'OuterPosition',[0  0  0.5  1/3])
plot(t/3600,rPwh,'r',t/3600,Pwh,'k','Linewidth',linewidth);
legend('set-point','P_{wh}',1,'Orientation','Vertical')
xlim([0 t_f/3600])
xlabel('Time [h]','FontSize',labelsize)
ylim([-1.5 1.5])

subplot(3,2,6)
set(gca,'FontSize',legsize,'OuterPosition',[0.5  0  0.5  1/3])
plot(t/3600,Prb,'k','Linewidth',linewidth);
legend('P_{rb}',1,'Orientation','Horizontal')
xlim([0 t_f/3600])
xlabel('Time [h]','FontSize',labelsize)
ylim([-2 2])