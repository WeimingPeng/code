clc
clear all
format long
warning off

controller='LQR';
estimator='UKF';
feedback='';
place = 'C:\Users\tereseva\Dropbox\Master\TereseCode\LQG\3percent\UKF\';
state = 'state';
input= 'input';
output2= 'output';

% **** System parameters ****
par = v4_new_4d_parameters();
% **************************
hours=4.5;
tfinal=hours*3600; % tfinal=6*3600;
t0=0;
Ts = 1;
t = t0:Ts:tfinal;
nf = length(t)-1;

load wG_wL_nm3percent.mat;

wG_in = 0.36*ones(1,nf);
wL_in = 8.64*ones(1,nf);
n_wG_in=n_wG_in(1:nf);
n_wL_in=n_wL_in(1:nf);
n_m=n_m(:,1:nf);

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

%##########################################
%%
x = zeros(4,nf+1);
x(:,1) = 0.001*x0;

y = zeros(3,nf);
y(:,1) = y0;
y_hat(:,1) = y0;
y_e = zeros(3,nf);
y_e(:,1) = 0;
P1_e = zeros(1,nf);
e = zeros(1,nf);
I = 0;
x_e = zeros(4,nf);
x_e(:,1) = 0*x0;

u_c = zeros(1,nf+1);
u_c2 = zeros(1,nf+1);

%% %%%UKF parameters%%%%
ff = 0.99;
L = 4;
alpha=1;
beta=0;
kapa=1;
lambda=alpha^2*(L+kapa)-L;
gama=sqrt(L+lambda);
Wc=[lambda/(L+lambda)+(1-alpha^2+beta) , ones(1,2*L)*lambda/(2*(L+lambda))];
Wm=[lambda/(L+lambda) , ones(1,2*L)*lambda/(2*(L+lambda))];
Wc1=ones(L,1)*Wc;
Wc2=ones(2,1)*Wc;
Wm1=ones(L,1)*Wm;
Wm2=ones(2,1)*Wm;

Pk_1 = (x_e(:,1) - x(:,1))*(x_e(:,1) - x(:,1))';
p11 = zeros(nf,L);
%% Intializibg the ode15s solver
MM=eye(4);
options=odeset('Mass',MM,'BDF','on','AbsTol',1e-10,'RelTol',1e-10,'MaxStep',1);


%***********************************************
%% The main loop
disp('********************************************')
disp('Simulation of the closed-loop system started')
disp('********************************************')

tic
for k=2:nf,
    if (k==2 || du_pc(k-1)~=0 || dwG_in(k-1)~=0 || dwL_in(k-1)~=0)
        disp(['New Operating point: Z=' num2str(u(1,k-1)) ' W_Gin=' num2str(u(2,k-1)) ' W_Lin=' num2str(u(3,k-1))])
        [x0,y0,par] =  v4_new_4d_initialize(u(:,k),par);
        [A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u(:,k-1),par);

        %#################################################
        F = expm(A*Ts);
        BL = F*(eye(size(A))-expm(-A*Ts))/A*B;
        H = C(sensors,:); %eq. 13.48 in Optimal State Estimation 
        M=eye(length(sensors)); %eq. 13.48 in Optimal State Estimation 
        DE = D(sensors,:); 

        B = BL(:,determ);
        LL = BL(:,stoch);

        DD = DE(:,determ);
        E = DE(:,stoch);

        % measurement noise vt := Ew+Mv
        % var (vt)=var(Ew+Mv)
        %         =E((Ew+Mv)(Ew+Mv)')
        %         =E(Ew+Mv)(w'E'+v'M')
        %         =E(Eww'E'+Ewv'M'+Mvw'E'+Mvv'M')
        %         =EQnE    +ENnM' +MNn'E'+ MRnM'
        %         =Eqn     +Hn  +Hn'  + MRnM'
        Hn = E * Nn * M';
        Nb = Qn * E' + Nn;
        R =  E * Qn * E' + Hn + Hn'+M*Rn*M';
        % Then incorporate LL
        Q = LL * Qn * LL';
        Nb = LL * Nb;

        % Enforce symmetry and check positivity
        Q_KF = (Q+Q')/2;
        R_KF = (R+R')/2;
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
    u_c(k-1)=-K_lqr*(x_e(:,k-1));   
    if(u_c(k-1)<-u_pc(k-1));u_c(k-1)=-u_pc(k-1);end
    if(u_c(k-1)>(1-u_pc(k-1)));u_c(k-1)=1-u_pc(k-1);end
    
    % =========== Simulation of System ===============
    u_in(:,k) = u_n(:,k-1)+[u_c(k-1);0;0];
    u_in_noNoise(:,k)=u(:,k-1)+[u_c(k-1);0;0];
    [tt,xt]=ode15s(@v4_new_4d_model,[k-1;k],x(:,k-1)+x0,options,u_in(:,k),'derivatives',par);
    x(:,k)= xt(end,:)'-x0;
    yt=v4_new_4d_model(k,xt(end,:)',u_in(:,k),'measurements',par);
    y(:,k) = yt(1:3)+n_m(:,k);
    %     [x(:,k),yt,k_]=ode45_fixed_step(FUN,t(k-1),Ts,x(:,k-1),u_in ,par,k_,odepar);
%     y(:,k) = yt+n_m(:,k);
     if mod(t(k),tfinal/10)==0 % Show message on every 10% of simulation
       disp([num2str(floor(t(k)/tfinal*100)) '% of simulation completed'])
    end
    
      
  % =========== Choosing Sigma Points ===============
    xhatk_1 = x_e(:,k-1)+x0;
    [RSK,p]=chol( (gama^2) * Pk_1 );
    if p==0,
        SK = RSK';
    else
             SK=sqrt(abs( (gama^2) * Pk_1 ));
    end

    Xk_1=[xhatk_1 , xhatk_1*ones(1,L)+SK , xhatk_1*ones(1,L)-SK ];

    %================== Propagation of Sigma Points =============
    Xk = zeros(L,2*L+1);
    for j=1:2*L+1,
        [tt,xt]=ode15s(@v4_new_4d_model,[k-1;k],Xk_1(:,j),options,u_in_noNoise(:,k),'derivatives',par);
        Xk(:,j) = xt(end,:)';    
    end

    xhatk_=sum((Xk).*Wm1,2);%/(2*L+1);

    temp1=Xk-xhatk_*ones(1,2*L+1);
    Pk_=Wc1.*temp1*temp1'+Q_KF;%
    %=============== Measurement Update ====================

    SK=chol( (gama^2) * Pk_ )';
    Xk=[xhatk_ , xhatk_*ones(1,L)+SK , xhatk_*ones(1,L)-SK];
    Yk = zeros(2,2*L+1);
    for j=1:2*L+1,
        y_k = v4_new_4d_model(t(k),Xk(:,j),u_in_noNoise(:,k),'measurements',par);
        Yk(:,j)= y_k(sensors);
    end
    
    yhatk_=sum(Wm2.*Yk,2);
    
    Pxkyk=Wc1.*(Xk-xhatk_*ones(1,2*L+1))*(Yk-yhatk_*ones(1,2*L+1))';%

    temp1=Yk-yhatk_*ones(1,2*L+1);
    Pyk_yk_=Wc2.*temp1*temp1'+R_KF;%


    Kk=Pxkyk/Pyk_yk_;

    xhatk=xhatk_+Kk*(y(sensors,k)-yhatk_);

    x_e(:,k)=xhatk-x0;
    y_hat(:,k)=v4_new_4d_model(t(k),xhatk,u_in_noNoise(:,k),'measurements',par);

    Pk=Pk_-Kk*Pyk_yk_*Kk';

    p11(k,:)=diag(Pk)';
    Pk_1=Pk/ff;

end
runningTime=toc;
