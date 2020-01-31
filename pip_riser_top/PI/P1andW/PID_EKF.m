clc
clear all
format long
warning off

controller='PI';
estimator='EKF';
feedback='';
place = 'C:\Users\tereseva\Dropbox\Master\TereseCode\PID\P1andW\EKF\';
state = 'state';
input= 'input';
output2= 'output';
% **** System parameters ****
par = v4_new_4d_parameters();
% **************************
hours=4.5;
tfinal=hours*3600;
t0=0;
Ts = 1;
t = t0:Ts:tfinal;
nf = length(t)-1;

load('wG_wL_nm.mat')

wG_in = 0.36*ones(1,nf);
wL_in = 8.64*ones(1,nf);

wG_in_n = wG_in + n_wG_in;
wL_in_n = wL_in + n_wL_in;

n_d = [n_wG_in; n_wL_in];

determ = 1; %deterministic input
stoch = [2 3]; %stochastic input
sensors = [1 3]; % measured outputs


%% --------------- Computing Expected Values for Kalman Filter --------
Qn = cov(n_d');
Rn = cov(n_m(sensors,:)');
N11 = cov(n_d(1,:),n_m(sensors(1),:)');
N12 = cov(n_d(1,:),n_m(sensors(2),:)');
N21 = cov(n_d(2,:),n_m(sensors(1),:)');
N22 = cov(n_d(2,:),n_m(sensors(2),:)');
Nn = [N11(1,2) N12(1,2) ; N21(1,2) N22(1,2) ];
%================================================

u_pc = 0.10*ones(1,nf);
u = [u_pc ; wG_in ; wL_in];
u_n = [u_pc ; wG_in_n ; wL_in_n];
disp('********************************')
disp('Initializing the state variables')
disp('********************************')
[x0,y0,par] =  v4_new_4d_initialize(u(:,1),par);
[A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u(:,1),par);
du_pc = [0 diff(u_pc)];
dwG_in = [0 diff(wG_in)];
dwL_in = [0 diff(wL_in)];

%############ PID Controller ##############
Kp = 1;
Ki = 0.0001;
%##########################################
%%
x = zeros(4,nf+1);
x(:,1) = 0.001*x0;
x_hat = zeros(4,nf);
x_hat(:,1) = 0.001*x0;
x_e=zeros(4,nf);
x_hat_pri=zeros(4,nf);

y = zeros(3,nf);
y(:,1) = y0;
y_hat = zeros(3,nf);
y_hat(:,1) = y0;
y_hat_pri=zeros(3,nf);
P1_hat = zeros(1,nf);
e = zeros(1,nf);
I = 0;
u_c = zeros(1,nf+1);
u_in=zeros(3,nf);
u_no_noise=zeros(3,nf);

Pk = (x_hat(:,1) - x(:,1))*(x_hat(:,1) - x(:,1))';
MM=eye(4);
options=odeset('Mass',MM,'BDF','on','AbsTol',1e-10,'RelTol',1e-10,'MaxStep',1);

%% The main loop
disp('********************************************')
disp('Simulation of the closed-loop system started')
disp('********************************************')
tic
for k=2:nf,

    if (k==2 || du_pc(k-1)~=0 || dwG_in(k-1)~=0 || dwL_in(k-1)~=0)
        disp(['New Operating point: Z=' num2str(u(1,k-1)) ' W_Gin=' num2str(u(2,k-1)) ' W_Lin=' num2str(u(3,k-1))])
        [x0,y0,par] =  v4_new_4d_initialize(u(:,k),par); %Linearize around working point, but not with control action.
        [A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u(:,k-1),par); % Continous system
        F = expm(A*Ts); %Discretize %eq, 4.17 in Linear System Theory and Design.
        BL = F*(eye(size(A))-expm(-A*Ts))/A*B;
        H = C(sensors,:);
        M = eye(length(sensors));
        DE = D(sensors,:);
        
        B = BL(:,determ);
        LL = BL(:,stoch);

        DD = DE(:,determ);
        E = DE(:,stoch);        

        % measurement noise vt := Ew+v
        Hn = E * Nn;
        Nb = Qn * E' + Nn;
        R = M*Rn*M'  + Hn + Hn' + E * Qn * E';
        % Then incorporate L
        Q = LL * Qn * LL';
        Nb = LL * Nb;

        % Enforce symmetry and check positivity
        Q = (Q+Q')/2;
        R = (R+R')/2;
        vr = real(eig(R));
        vqnr = real(eig([Q Nb;Nb' R]));
        if min(vr)<0 || (Ts==0 && min(vr)==0),
           error(...
             ['The covariance matrix\n' ...
              '    E{(Hw+v)(Hw+v)''} = [H,I]*[Qn Nn;Nn'' Rn]*[H'';I]\n'...
              'must be positive definite.'])
        elseif min(vqnr)<-1e2*eps*max(0,max(vqnr)),
            warning(1,'The matrix [G 0;H I]*[Qn Nn;Nn'' Rn]*[G 0;H I]'' should be positive semi-definite.')
        end

    end
    
    % ============== PID Controller ==================
    P1_hat(k-1)=y_hat(1,k-1);
    e(k-1) = P1_hat(k-1) - y0(1);
    I = I + e(k-1)*Ts;
    u_c(k-1) = Kp*(e(k-1)+Ki*I);    
    if(u_c(k-1)<-u_pc(k-1));u_c(k-1)=-u_pc(k-1);end
    if(u_c(k-1)>(1-u_pc(k-1)));u_c(k-1)=1-u_pc(k-1);end
    
     % =========== Simulation of System ===============
     
     u_in(:,k) = u_n(:,k-1)+[u_c(k-1);0;0];
     u_no_noise(:,k)=u(:,k-1)+[u_c(k-1);0;0];

    [tt,xt]=ode15s(@v4_new_4d_model,[k-1;k],x(:,k-1)+x0,options,u_in(:,k),'derivatives',par);
    x(:,k)= xt(end,:)'-x0;
    yt=v4_new_4d_model(k,xt(end,:)',u_in(:,k),'measurements',par);
    y(:,k) = yt(1:3)+n_m(:,k);

    if mod(t(k),tfinal/10)==0 % Show message on every 10% of simulation
       disp([num2str(floor(t(k)/tfinal*100)) '% of simulation completed'])
    end
    
    % =========== Extended Kalman Filter ===============
    % p.409 in Optimal State Estimation. By Dan J. Simon

    Pk_pri = F*Pk*F' + Q;
    
    [tt,xt]=ode15s(@v4_new_4d_model,[k-1;k],x_hat(:,k-1)+x0,options,u_no_noise(:,k),'derivatives',par);
    x_hat_pri(:,k)= xt(end,:)'-x0;
    yt(1:3,k)=v4_new_4d_model(k,xt(end,:)',u_no_noise(:,k),'measurements',par);
    y_hat_pri(:,k) = yt(1:3,k);
        
    K = Pk_pri*H'/(H*Pk_pri*H'+R);
    x_hat(:,k)= x_hat_pri(:,k) + K*(y(sensors,k) - y_hat_pri(sensors,k));
    x_e(:,k)=x_hat(:,k)-x0;
    yt(1:3,k) = v4_new_4d_model(t(k),x_hat(:,k)+x0,u_no_noise(:,k),'measurements',par);
    y_hat(:,k) = yt(1:3,k);
    Pk = (eye(size(Pk_pri))-K*H)*Pk_pri;

end
runningTime=toc