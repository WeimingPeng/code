clc
clear
% **** System parameters ****
par = v4_new_4d_parameters();
% **** System Inputs ****
z0 = 0.05;
wG_in0 = 0.3;
wL_in0 = 5.64;
global u0;
u0 = [z0;wG_in0;wL_in0];
% ******** Finding initial states *************
disp('Initializing...')
tic
[x0,y0,par] =  v4_new_4d_initialize(u0,par);
toc
disp('Initialization completed.');
disp('#########################');
%% ******** Linearization *************
% options = [0,1e-5,1e-5,1e-12,0,0,0,0,0,0,0,0,0,1e5];
% [x0,u0,y,dx,options] = trim('v2_New_4D_mdl',x0,u0,y0,[],[],[],[],[],options);
para = [1e-5 0 0];
[A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u0,par);
disp('Linearized model in operating point Z=0.05')
sys = ss(A,B,C,D);
set(sys,'inputname',{'Z','wG_in','wL_in'},'outputname',{'P1','P2','W'},'statename',{'m_G1','m_L1','m_G2','m_L2'});
present(sys)
G1 = tf(sys('P1','Z'));
G2 = tf(sys('P2','Z'));
G3 = tf(sys('W','Z'));
present(zpk(G1))
present(zpk(G2))
present(zpk(G3))

G10 = abs(freqresp(G1,0))
G20 = abs(freqresp(G2,0))
G30 = abs(freqresp(G3,0))

disp('poles:')
p = pole(G1);
disp(p)
disp('Period of the linearized model [sec]')
disp(2*pi/abs(imag(p(3))))

disp('zeros of P1:')
disp(zero(G1))
disp('zeros of P2:')
disp(zero(G2))
disp('zeros of W:')
disp(zero(G3))

%% **********************
tic
disp('Simulation of nonlinear model');
%dt = 0.005;    % Simulation step time (sec)
t_f = 1;        % Final simulation time (hour)

dx = [-3;0;0;0];
u=[z0;wG_in0;wL_in0];
M=eye(4);
options=odeset('Mass',M,'BDF','on','AbsTol',1e-12,'RelTol',1e-12,'MaxStep',1);
[ts,xs]=ode15s(@v4_new_4d_model,[0;3600*t_f],x0 + dx,options,u0,'derivatives',par);
% *** Calculating measurements ****
n = length(ts);
ys = zeros(3,n);
for k = 1:n
    yt = v4_new_4d_model(ts(k),xs(k,:),u0,'measurements',par);
    ys(:,k) = yt;
end

disp('Simulation completed.');
toc
disp('#########################');
%% ******** Simulation Using ode45 fixed-steo ********
% disp('Simulation using ode45 fixed step');
% tic
% Ts = 0.05;
% tfinal = t_f*3600;
% nf = t_f*3600/Ts;
% t = 0:Ts:t_f*3600;
% L = length(x0);
% x = zeros(L,nf+1);
% x(:,1) = x0+dx;
% y = zeros(3,nf);
% y(:,1) = v4_new_4d_model(t(1),x(:,1),u0,'measurements',par);
% % Intializibg the ode solver
% odepar = ode45('Fehlberg',1e-3,Ts,1,0,0);
% %Fehlberg
% %Dormand-Prince
% initialize_solver;
% if (odepar.pair==0)
%     k_sp = zeros(L,7);
% elseif(odepar.pair==1)
%     k_sp = zeros(L,6);
% end
% 
% for k=2:nf
%     [x(:,k),y(:,k),k_]=ode45_fixed_step(FUN,t(k-1),Ts,x(:,k-1),u0 ,par,k_,odepar);
%     if mod(t(k-1),tfinal/10)==0 % Show message on every 10% of simulation
%        disp([num2str(floor(t(k-1)/tfinal*100)) '% of simulation completed'])
%     end
% end
% disp('Simulation completed.');
% toc
% disp('#########################');
% %% ********************* Figure 1 **************************
% figure(1)
% clf
% set(gcf,'Color',[1,1,1])
% set(gca,'FontSize',8)
% 
% subplot(3,1,1)
% set(gca,'FontSize',8,'OuterPosition',[0  2/3  1  1/3])
% plot(ts/60,ys(1,:),'k');
% hold on
% plot(t(1:k)/60,y(1,:),'--r');
% title('a) Pressure at inlet of pipeline (\itP\rm\bf_1)','FontSize',10,'FontWeight','bold','FontName','Times');
% ylabel('\itP\rm_1 [bar]','FontSize',10,'FontName','Times')
% xlabel('time [min]','FontSize',10,'FontName','Times')
% xlim([0 t_f*60])
% %ylim([68 80])
% 
% 
% subplot(3,1,2)
% set(gca,'FontSize',8,'OuterPosition',[0  1/3  1  1/3])
% plot(ts/60,ys(2,:),'k');
% hold on
% plot(t(1:k)/60,y(2,:),'--r');
% title('b) Pressure at top of riser (\itP\rm\bf_2)','FontSize',10,'FontWeight','bold','FontName','Times');
% ylabel('\itP\rm_2 [bar]','FontSize',10,'FontName','Times')
% xlabel('time [min]','FontSize',10,'FontName','Times')
% xlim([0 t_f*60])
% %ylim([50 62])
% 
% 
% subplot(3,1,3)
% set(gca,'FontSize',8,'OuterPosition',[0  0  1  1/3])
% plot(ts/60,ys(3,:),'k');
% hold on
% plot(t(1:k)/60,y(3,:),'--r');
% title('c) Mass flow rate form choke valve (\itw_{out}\rm\bf)','FontSize',10,'FontWeight','bold','FontName','Times');
% ylabel('\itw_{out}\rm [kg/s]','FontSize',10,'FontName','Times')
% xlabel('time [min]','FontSize',10,'FontName','Times')
% xlim([0 t_f*60])
% %ylim([5 35])
