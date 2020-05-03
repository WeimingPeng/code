% Simulation of the new 4d model
clc
clear all
% **** System parameters ****
par = v2_new_5d_parameters();
% **** System Inputs ****
z1 = 30;
z2 = 0.04;
global u0;
u0 = [z1;z2;0;0];
% ******** Finding initial states *************
disp('Initializing...')
tic
[x0,y0,par] =  v2_new_5d_initialize(u0,par);
x0
toc
disp('Initialization completed.');
disp('#########################');
para = [1e-6;0;1];
%[A,B,C,D]=linmod('linerazation_simulink',x0,u0,para);
[A,B,C,D] = v2_new_5d_linmodl_num(x0,u0,par);
sys = ss(A,B,C,D);
pole(sys)
%dt = 0.005;    % Simulation step time (sec)
t_f = 2;        % Final simulation time (hour)

%% **********************
disp('Simulation of nonlinear model');
dx = [0;-0.02;0;0;0];
u=[z1;z2];
M=eye(5);
options=odeset('Mass',M,'BDF','on','AbsTol',1e-10,'RelTol',1e-10,'MaxStep',1);
[t,x]=ode15s(@v2_new_5d_model,[0;3600*t_f],x0 + dx,options,u0,'derivatives',par);
% *** Calculating measurements ****
n = length(t);
y = zeros(6,n);
for k = 1:n
    yt = v2_new_5d_model(t(k),x(k,:),u0,'measurements',par);
    y(:,k) = yt([1 2 3 4 7 9]);
end

disp('Simulation completed.');
disp('#########################');

%% ********************* Figure 1 **************************
figure(1)
clf
set(gcf,'Color',[1,1,1])
set(gca,'FontSize',8)

subplot(3,2,1)
set(gca,'FontSize',8,'OuterPosition',[0  2/3  1/2  1/3])
plot(t/60,y(1,:),'k');
title('a) Bottom-hole pressure (\itP\rm\bf_{bh})','FontSize',12,'FontWeight','bold','FontName','Times');
ylabel('\itP\rm_{bh} [bar]','FontSize',12,'FontName','Times')
xlabel('time [min]','FontSize',12,'FontName','Times')
xlim([0 t_f*60])
%ylim([68 80])

subplot(3,2,2)
set(gca,'FontSize',8,'OuterPosition',[1/2  2/3  1/2  1/3])
plot(t/60,y(2,:),'k');
title('b) Well-head pressure (\itP\rm\bf_{wh})','FontSize',12,'FontWeight','bold','FontName','Times');
ylabel('\itP\rm_{wh} [bar]','FontSize',12,'FontName','Times')
xlabel('time [min]','FontSize',12,'FontName','Times')
xlim([0 t_f*60])
%ylim([5 35])

subplot(3,2,3)
set(gca,'FontSize',8,'OuterPosition',[0  1/3  1/2  1/3])
plot(t/60,y(3,:),'k');
title('c) Inlet mass flow rate (\itw\rm\bf_{in})','FontSize',12,'FontWeight','bold','FontName','Times');
ylabel('\itw\rm_{in} [kg/s]','FontSize',12,'FontName','Times')
xlabel('time [min]','FontSize',12,'FontName','Times')
xlim([0 t_f*60])
%ylim([50 62])

subplot(3,2,4)
set(gca,'FontSize',8,'OuterPosition',[1/2  1/3  1/2  1/3])
plot(t/60,y(4,:),'k');
title('d) Pressure at inlet of pipeline (\itP\rm\bf_1)','FontSize',12,'FontWeight','bold','FontName','Times');
ylabel('\itP\rm_1 [bar]','FontSize',12,'FontName','Times')
xlabel('time [min]','FontSize',12,'FontName','Times')
xlim([0 t_f*60])
%ylim([50 62])

subplot(3,2,5)
set(gca,'FontSize',8,'OuterPosition',[0  0 1/2  1/3])
plot(t/60,y(5,:),'k');
title('e) Pressure at top of riser (\itP\rm\bf_2)','FontSize',12,'FontWeight','bold','FontName','Times');
ylabel('\itP\rm_2 [bar]','FontSize',12,'FontName','Times')
xlabel('time [min]','FontSize',12,'FontName','Times')
xlim([0 t_f*60])
%ylim([5 35])


subplot(3,2,6)
set(gca,'FontSize',8,'OuterPosition',[1/2  0  1/2  1/3])
plot(t/60,y(6,:),'k');
title('f) Outlet mass flow rate (\itw\rm\bf_{out})','FontSize',12,'FontWeight','bold','FontName','Times');
ylabel('\itw\rm_{out} [kg/s]','FontSize',12,'FontName','Times')
xlabel('time [min]','FontSize',12,'FontName','Times')
xlim([0 t_f*60])
%ylim([5 35])
disp([min(y(1,:)) max(y(1,:))])
disp([min(y(2,:)) max(y(2,:))])
disp([min(y(3,:)) max(y(3,:))])
disp([min(y(4,:)) max(y(4,:))])
disp([min(y(5,:)) max(y(5,:))])
disp([min(y(6,:)) max(y(6,:))])
