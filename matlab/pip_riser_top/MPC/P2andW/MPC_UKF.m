clc
clear all
format long


controller='MPC';
estimator='UKF';
feedback='ff098';
place = 'C:\Users\tereseva\Dropbox\Master\TereseCode\MPC\3percent\UKF\';
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


load('wG_wL_nm3percent.mat')

n_wG_in = n_wG_in(1:nf);
n_wL_in = n_wL_in(1:nf);
n_m=n_m(:,1:nf);

wG_in = 0.36*ones(1,nf);
wL_in = 8.64*ones(1,nf);

wG_in_n = wG_in + n_wG_in;
wL_in_n = wL_in + n_wL_in;

n_d = [n_wG_in; n_wL_in];

determ = 1; %deterministic
stoch = [2 3];%stochastic
sensors =[2 3];

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
[A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u(:,1),par);
du_pc = [0 diff(u_pc)];
dwG_in = [0 diff(wG_in)];
dwL_in = [0 diff(wL_in)];


%% Initialize MPC

%%Dimensions
Hu=15;  %Control horizon
Hp=Hu;  %Prediction horizon
Hw=1;   %When you start controlling


%Table 3.1 in Predictive Control with constraints by Maciejowski
l=1; %inputs
n=4; %States
m=4; %Controlled outputs ,z=cx, c=eye(4);
Q_MPC=zeros(m*(Hp-Hw+1),m*(Hp-Hw+1));
R_MPC=zeros(l*Hu,l*Hu);
Psi=zeros(m*(Hp-Hw+1),n);
Gamma=zeros(m*(Hp-Hw+1),l);
Theta=zeros(m*(Hp-Hw+1),l*Hu);


%Initialize constraints matrix
Ff   = zeros(2*Hu,Hu);
%f1   = zeros(2*Hu,1);
f    = zeros(2*Hu,1);
Gcons = zeros(2*Hu,Hu);
g1    = zeros(2*Hu,1);
%g     = zeros(2*Hu,1);
W    = zeros(2*Hu,Hu);
w    = zeros(2*Hu,1);

%Weighting on the different states and input.
Q_MPCd=1*diag([1 0 1 0]);
R_MPCd=1e5;
Q_MPC=blkdiag2(Q_MPCd,Hu);
R_MPC=blkdiag2(R_MPCd,Hu);



%% Constraints

bias=u(1,1);
uMax=1;  %uMin<u<uMax
uMin=0; 

uMax=uMax-bias;
uMin=uMin-bias;

violation=0.5;
zMax1=violation*x0(1);
zMin1=-violation*x0(1);
zMax2=violation*x0(2);
zMin2=-violation*x0(2);
zMax3=violation*x0(3);
zMin3=-violation*x0(3);
zMax4=violation*x0(4);
zMin4=-violation*x0(4);
uRateMax=1e-1;
uRateMin=-1e-1;


%% Ff is 1. line on left side of constraint (P)
for j = 1:Hu,
    f(2*(j-1)+1,1) = -uMax;
    f(2*j,1)       = uMin;    
    for i = j:Hu,
        Ff(2*i-1,j)       = 1;
        Ff(2*i,(j-1)+1:j) = -1;
    end
end

f1    = Ff(1:2*Hu,1);

constZ=[-zMax1; zMin1;-zMax2;zMin2;-zMax3;zMin3;-zMax4;zMin4];
g=repmat(constZ,Hu,1);

for i = 1:Hu*m,
    Gcons(2*i-1,i) = 1;
    Gcons(2*i,i)   = -1;
end

%% W and w matrices in constraints
%% W is 3. line on left side of constraint
for i = 1:Hu,
    w(2*i-1,1) = uRateMax;
    w(2*i,1)   = -uRateMin;
    W(2*i-1,i) = 1;
    W(2*i,i)   = -1;
end
   
x = zeros(4,nf+1);
x(:,1) = 0.001*x0;
y = zeros(size(y0,1),nf);
y(:,1) = y0;
y_hat = zeros(3,nf);
y_hat(:,1)= y0;
x_e = zeros(4,nf);
x_e(:,1) = 0*x0;

u_in=zeros(3,nf);
u_in_noNoise = zeros(3,nf);
u_MPC=zeros(1,nf);

%% %%%UKF parameters%%%%
ff = 0.98;
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
p11 = zeros(nf,4);
I=0;
MM=eye(length(sensors));
M=eye(4);
%% Initializing Ode15s solver and QPsolver
options=odeset('Mass',M,'BDF','on','AbsTol',1e-10,'RelTol',1e-10,'MaxStep',1);
optionsQP = optimset('Display','off','LargeScale','off','Algorithm','active-set');
%***********************************************
%% The main loop
disp('********************************************')
disp('Simulation of the closed-loop system started')
disp('********************************************')

tic
for k=2:nf

    if (k==2 || du_pc(k-1)~=0 || dwG_in(k-1)~=0 || dwL_in(k-1)~=0)
        disp(['New Operating point: Z=' num2str(u(1,k-1)) ' W_Gin=' num2str(u(2,k-1)) ' W_Lin=' num2str(u(3,k-1))])
        %[x0,y0,par] =  v4_new_4d_initialize(u(:,k),par); %Linearize around working point, but not with control action.
        [Ac,Bc,Cc,Dc]= v4_new_4d_linmod_knut(x0,y0,u(:,k-1),par); %Index c: continiuous
        %[Ac,Bc,Cc,Dc] = v3_new_4d_LinModel_num(x0,y0,u(:,k-1),par);
        sysc = ss(Ac,Bc,Cc(1:3,:),Dc(1:3,:));
        sysd = c2d(sysc,Ts);
        A = sysd.a; B = sysd.b; C = sysd.c; D = sysd.d;

        BB = B(:,determ);
        Bd= B(:,stoch);
        C = C(sensors,:); %eq. 13.48 in Optimal State Estimation 
        D = D(sensors,:); 
        
        DD=D(:,determ);
        E=D(:,stoch);
        
        Hn=E*Nn*MM';
        Nb=Qn*E'+Nn;
        R_KF=E*Qn*E'+Hn+Hn'+MM*Rn*MM';
        
        Q_KF=Bd*Qn*Bd';
        Nb=Bd*Nb;
        
        % Enforce symmetry and check positivity
        Q_KF = (Q_KF+Q_KF')/2;
        R_KF = (R_KF+R_KF')/2;
        vr = real(eig(R_KF));
        vqnr = real(eig([Q_KF Nb;Nb' R_KF]));
        if min(vr)<0 || (Ts==0 && min(vr)==0),
           error(...
             ['The covariance matrix\n' ...
              '    E{(Hw+v)(Hw+v)''} = [H,I]*[Qn Nn;Nn'' Rn]*[H'';I]\n'...
              'must be positive definite.'])
        elseif min(vqnr)<-1e2*eps*max(0,max(vqnr)),
            warning(1,'The matrix [G 0;H I]*[Qn Nn;Nn'' Rn]*[G 0;H I]'' should be positive semi-definite.')
        end
        
    
    %% =========== MPC ===============

        an=A;
        %Make Psi and Gamma matrices
        Gamma(1:n,1:l)=BB;
        GammaDist(1:n,1:2)=Bd;
        for i=1:n:Hp*n
            Psi(i:i+(n-1),1:n)=an;            
            if i> 1
                Gamma(i:(i+(n-1)),1:l)=Gamma(i-n:(i-1),1:l)+Psi(i-n:(i-1),1:n)*BB;
            end
            an=an*A;
        end
        
        %Dummy matrix to construct Theta.        
        matrix(1:n,1:n)=eye(size(A));
        for i=n+1:n:Hu*n            
            matrix(i:i+n-1,1:n)=A*matrix(i-n:i-1,1:n)+eye(m);
        end     

        matrixB(1:n,1:l)=matrix(1:n,1:n)*BB;
        for i=1:n:Hu*n
            matrixB(i:i+n-1,1:l)=matrix(i:i+n-1,1:4)*BB;
        end

        % Make Theta matrix
        Theta(1:n,:)=[matrixB(1:n,1:l) zeros(n,Hu-l)];
        for i=5:n:Hu*n
              Theta(i:i+n-1,:)=[matrixB(i:i+n-1,1:l) Theta(i-n:i-1,1:Hu-l)];
        end
     end   
    %% Formulate H and G.
    x_prev=x_e(:,k-1);
    u_control_prev=u_MPC(k-1);

    T=repmat([0 0 0 0]',Hu,1);
    EPS=T-Psi*x_prev-Gamma*u_control_prev; %eq. 3.6      
    H=Theta'*Q_MPC*Theta+R_MPC;     %eq. 3.12 
    H=(H+H')/2;                     %To ensure symmetry/compensate for numerical inaccuracies
    G=2*Theta'*Q_MPC*EPS;           %eq. 3.11



    %% Formulate constraints as eq. 3.41

%         Ff                %1st line, left side   
     GcFi=Gcons*Theta; %2nd line, left side 
%         W                 %3rd line, left side

    %Eps=-Psi*x_prev-Gamma*u_control_prev; %2nd line right side, midle part 

    ohm1=-f1*u_control_prev-f;      %1st line, right side
    g1=Gcons*EPS-g;               %2nd line, right side
  %  w                             %3rd line, right side

    OhmL=[Ff;GcFi;W];      % Total vector, left side of constraints
    OhmR = [ohm1; g1; w];  % Total vector, right side of constraints
                           % OhmL*deltaU<=OhmR, evt. Ax<=b  
        
     %% Quadratic program, constrained case:
     
    [delU(:,k),fval,exFlag,output,lambda] = quadprog(H, -0.5*G, OhmL, OhmR, [],[],[],[],[],optionsQP);
    
     
    %% Input to simulation of plant

    u_MPC(k) = delU(1,k)+u_MPC(k-1); % Add up to total control (deltaU computed)
    if (u_MPC(k)>1-u(1,k))
        u_MPC(k)=1-u(1,k);
    elseif(u_MPC(k)<-u(1,k))
        u_MPC(k)=-u(1,k);
    end
    
    u_in(:,k)=[u_n(1,k)+u_MPC(k); u_n(2,k); u_n(3,k)];  
    u_in_noNoise(:,k)=[u(1,k)+u_MPC(k); u(2,k); u(3,k)];
    
%% =========== Simulation ===============
    [tt,xt]=ode15s(@v4_new_4d_model,[k-1;k],x(:,k-1)+x0,options,u_in(:,k),'derivatives',par);
    x(:,k) = xt(end,:)'-x0;
    yt = v4_new_4d_model(t(k),x(:,k)+x0,u_in(:,k),'measurements',par);
    y(:,k) = yt(1:3)+ n_m(:,k);
   
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

    xhatk_=sum((Xk).*Wm1,2);

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
    Pyk_yk_=Wc2.*temp1*temp1'+R_KF;


    Kk=Pxkyk/Pyk_yk_;

    xhatk=xhatk_+Kk*(y(sensors,k)-yhatk_);

    x_e(:,k)=xhatk-x0;
    y_hat(:,k)=v4_new_4d_model(t(k),xhatk,u_in_noNoise(:,k),'measurements',par);

    Pk=Pk_-Kk*Pyk_yk_*Kk';

    p11(k,:)=diag(Pk)';
    Pk_1=Pk/ff;
  
end
%%
runningTime=toc;