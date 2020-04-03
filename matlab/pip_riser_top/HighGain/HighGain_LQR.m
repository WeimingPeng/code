clc
clear

controller='LQG';
estimator='HighGain';
feedback='';
place = 'C:\Users\tereseva\Dropbox\Master\TereseCode\NewHighGain\';
state = 'state';
input= 'input';
output2= 'output';

% **** System parameters ****
par = v4_new_4d_parameters();
% **** System Inputs ****
z0 = 0.1;
wG_in0 = 0.36;
wL_in0 = 8.64;
u0 = [z0;wG_in0;wL_in0];
%% ******** Linearization *************
t_f = 4.5; %hours
tfinal=t_f*3600; %seconds
t0=0;
Ts = 1; % step
t = t0:Ts:tfinal;
nf = length(t)-1;

load wG_wL_nm.mat

n_wG_in = n_wG_in(1:nf);
n_wL_in = n_wL_in(1:nf);

n_m = n_m(:,1:nf);

n_d = [n_wG_in; n_wL_in];

determ = 1; %deterministic input
stoch = [2 3]; %stochastic input
sensors = 2; % measured outputs


u1 = z0*ones(1,nf);
u2 = wG_in0*ones(1,nf);
u3 = wL_in0 *ones(1,nf);

u = [u1 ; u2 ; u3];
u_n = [u1 ; u2+n_wG_in; u3+n_wL_in];
disp('********************************')
disp('Initializing the state variables')
disp('********************************')
tic
[x0,y0,par] =  v4_new_4d_initialize(u0,par);
toc
[A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u0,par);
du1 = [0 diff(u1)];

%##########################################
%%
x = zeros(4,nf+1);
dx = [-3;0;0;0];
x(:,1) = x0+dx;

y = zeros(3,nf);
y(:,1) = y0(1:3);
y_e = zeros(3,nf);
y_e(:,1) = y(:,1);

x_e = zeros(4,nf);
x_e(:,1) = x0-dx;

z_e = zeros(4,nf);
z_e(:,1) = [x_e(1,1);x_e(2,1);1e5*y_e(3,1);x_e(4,1)];
%############ Controller Weight Matrices ##############
Q1 = 1e-3*diag([1 0 1 0]);
R1 = 1e3;
N1 = zeros(6,1);
x_r = x0;
u_c = zeros(1,nf);
%##########################################
%%%%% High-Gain Observer Parameter %%%%
ep = 5e-6;
%% Intializibg the ode solver
M=eye(4);
options=odeset('Mass',M,'BDF','on','AbsTol',1e-10,'RelTol',1e-10,'MaxStep',1);
%***********************************************
%% The main loop
disp('********************************')
disp('Simulation of the system started')
disp('********************************')
tic
for k=2:nf
    if (k==2 || du1(k-1)~=0)
        disp(['New Operating point: Z=' num2str(u(1,k-1)) ' W_Gin=' num2str(u(2,k-1)) ' W_Lin=' num2str(u(3,k-1))])
        %[x0,y0,par] =  v1_new_6d_initialize(u(:,k),par);
        [A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u0,par);
        
        Ob = obsv(A,C(sensors,:));
        Co = ctrb(A,B(:,[1 2]));
        
        % Number of unobservable states
        unob = rank(A)-rank(Ob);
        if unob>0
            disp('system is not observable');
        else
            disp('system is observable');
        end
        % Number of uncontrollable states
        unco=rank(A)-rank(Co);
        if unco>0
            disp('system is not controllable');
        else
            disp('system is controllable');
        end
        %#################################################
        F = expm(A*Ts); %Discretize %eq. 4.17 in Linear System Theory and Design.
        BL  = F*(eye(size(A))-expm(-A*Ts))/A*B; %inv(A)(F-I)B
        H = C(sensors,:); %eq. 13.48 in Optimal State Estimation 
        M=eye(length(sensors)); %eq. 13.48 in Optimal State Estimation 
        DE = D(sensors,:); 

        B = BL(:,determ);
        LL = BL(:,stoch);

        DD = DE(:,determ);
        E = DE(:,stoch);


        [S,E,K] = dare(F,B(:,1),Q1,R1);
    end
    
    % ============== LQR Controller ==================;   
    u_c(k-1) = -(k>0)*K*(x_e(:,k-1)-x_r);   
    if(u_c(k-1)<-u1(k-1));u_c(k-1)=-u1(k-1);end
    if((u_c(k-1)+u1(k-1))>1);u_c(k-1)=1-u1(k-1);end
    % =========== Simulation of System ===============
    u_in(:,k) = u_n(:,k-1)+[u_c(k-1);0;0];
    [~,xt]=ode15s(@v4_new_4d_model,[k-1;k],x(:,k-1),options,u_in(:,k),'derivatives',par);
    x(:,k) = xt(end,:)';
    yt = v4_new_4d_model(t(k),x(:,k),u_in(:,k),'measurements',par);
    y(:,k) = yt(1:3) + n_m(:,k);
    if mod(t(k),tfinal/10)==0 % Show message on every 10% of simulation
        disp([num2str(floor(t(k)/tfinal*100)) '% of simulation completed'])
    end
    
    %================== High-Gain Observer =============
    u_in(:,k) = u(:,k-1)+[u_c(k-1);0;0]; 
    [~,zt]=ode15s(@v1_new_4d_Observer,[k-1;k],z_e(:,k-1),options,u_in(:,k),y(sensors,k),'derivatives',par,ep);
    z_e(:,k) = zt(end,:)';
    x_e(:,k) = v1_new_4d_Observer(t,z_e(:,k),u_in(:,k),y(sensors,k),'measurements',par,ep);

    y_k = v4_new_4d_model(t(k),x_e(:,k),u_in(:,k),'measurements',par);
    y_e(:,k)=y_k(1:3);
    
        
end
runningTime=toc;
y_hat=y_e;
