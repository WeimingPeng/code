clc
clear all
format long
warning off

controller='LQG';
estimator='EKF';
feedback='';
place = 'C:\Users\tereseva\Dropbox\Master\TereseCode\LQG\3percent\EKF\';
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

load wG_wL_nm3percent.mat;

wG_in = 0.36*ones(1,nf);
wL_in = 8.64*ones(1,nf);

wG_in_n = wG_in + n_wG_in;
wL_in_n = wL_in + n_wL_in;

n_d = [n_wG_in; n_wL_in];

determ = 1; %deterministic input
stoch = [2 3]; %stochastic input
sensors = [2 3]; % measured outputs
   
%% --------------- Computing Expected Values for Kalman Filter --------
Qn = cov(n_d');
Rn = cov(n_m(sensors,:)');
nn11 = cov(n_d(1,:),n_m(sensors(1),:)');
nn12 = cov(n_d(1,:),n_m(sensors(2),:)');
nn21 = cov(n_d(2,:),n_m(sensors(1),:)');
nn22 = cov(n_d(2,:),n_m(sensors(2),:)');
Nn = [nn11(1,2) nn12(1,2) ; nn21(1,2) nn22(1,2) ];
%================================================

u_pc = 0.10*ones(1,nf);
u_pc = u_pc(1:nf);
u = [u_pc ; wG_in ; wL_in];
u_n = [u_pc ; wG_in_n ; wL_in_n];
disp('********************************')
disp('Initializing the state variables')
disp('********************************')
[x0,y0,par] =  v4_new_4d_initialize(u(:,1),par);

%%
x = zeros(4,nf+1);
x(:,1) = 0.001*x0;
x_hat = zeros(4,nf);
x_hat(:,1) = 0.001*x0;


x_hat_pri=zeros(4,nf);
y_hat_pri=zeros(3,nf);

y = zeros(3,nf);
y(:,1) = y0;
y_hat = zeros(3,nf);
y_hat(:,1) = v4_new_4d_model(0,x0,u(:,1),'measurements',par);

% Initialize 
%***********************************************
%% The main loop
disp('********************************************')
disp('Simulation of the closed-loop system started')
disp('********************************************')
P = (x_hat(:,1) - x(:,1))*(x_hat(:,1) - x(:,1))';
p11 = zeros(nf,4);
u_c = zeros(3,nf+1);

MM=eye(4);
% Initialize solver
options=odeset('Mass',MM,'BDF','on','AbsTol',1e-10,'RelTol',1e-10,'MaxStep',1);
tic
for k=2:nf,

    if (k==2 || du_pc(k-1)~=0 || dwG_in(k-1)~=0 || dwL_in(k-1)~=0)
        disp(['New Operating point: Z=' num2str(u(1,k-1)) ' W_Gin=' num2str(u(2,k-1)) ' W_Lin=' num2str(u(3,k-1))])
        [x0,y0,par] =  v4_new_4d_initialize(u(:,k),par);
        [A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u(:,k-1),par);
        F = expm(A*Ts); %Discretize %eq. 4.17 in Linear System Theory and Design.
        BL  = F*(eye(size(A))-expm(-A*Ts))/A*B; %inv(A)(F-I)B
        H = C(sensors,:); %eq. 13.48 in Optimal State Estimation 
        M=eye(length(sensors)); %eq. 13.48 in Optimal State Estimation 
        DE = D(sensors,:); 

        B = BL(:,determ);
        LL = BL(:,stoch);

        DD = DE(:,determ);
        E = DE(:,stoch);
        
        Hn = E * Nn * M';
        Nb = Qn * E' + Nn;
        R =  E * Qn * E' + Hn + Hn'+M*Rn*M';
        % Then incorporate LL
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
        
        

 % ============== LQR Controller ==================
        Q_lqr=1e-3*diag([1 0 1 0]);
        R_lqr=1e9;
        N_lqr=zeros(4,1);    
        [K_lqr,S,e] =dlqr(F,B,Q_lqr,R_lqr,N_lqr);
    end
    
    u_c(k-1)=-K_lqr*(x_hat(:,k-1));
    if(u_c(k-1)<-u_pc(k-1));u_c(k-1)=-u_pc(k-1);end
    if(u_c(k-1)>(1-u_pc(k-1)));u_c(k-1)=1-u_pc(k-1);end
    
    % =========== Simulation of System ===============
     
     
    u_in(:,k)=[u(1,k)+u_c(k-1); u_n(2,k); u_n(3,k)];  
    u_no_noise(:,k)=[u(1,k)+u_c(k-1); u(2,k); u(3,k)];
     
     [tt,xt]=ode15s(@v4_new_4d_model,[k-1;k],x(:,k-1)+x0,options,u_in(:,k),'derivatives',par);
    x(:,k)= xt(end,:)'-x0;
    yt=v4_new_4d_model(k,xt(end,:)',u_in(:,k),'measurements',par);
    y(:,k) = yt(1:3)+n_m(:,k);

    
    % =========== Extended Kalman Filter ===============
    % p.409 in Optimal State Estimation. By Dan J. Simon
    %3a) 
    %F =expm(A*Ts)
    %BL=(F-I)inv(A)B
    %L =BL(:,sensors)
  
    %b)
    P_pri=F*P*F'+Q;
     
    [tt,xt]=ode15s(@v4_new_4d_model,[k-1;k],x_hat(:,k-1)+x0,options,u_no_noise(:,k),'derivatives',par);
    x_hat_pri(:,k)= xt(end,:)'-x0;
    yt=v4_new_4d_model(k,xt(end,:)',u_no_noise(:,k),'measurements',par);
    y_hat_pri(:,k) = yt(1:3);
    % y_hat_pri is calculated based on x_hat_pri.
    
    %c)
    %H=C(sensors,:);
    %M=eye(length(sensors));
    
    %d) 
    K=P_pri*H'/(H*P_pri*H'+R);
    x_hat(:,k)= x_hat_pri(:,k) + K*(y(sensors,k) - y_hat_pri(sensors,k));
    y_hat(:,k) = v4_new_4d_model(t(k),x_hat(:,k)+x0,u_no_noise(:,k),'measurements',par);
    P=(eye(size(P_pri))-K*H)*P_pri;
 
end
runningTime=toc;