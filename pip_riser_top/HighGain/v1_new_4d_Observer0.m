function [zdot,xhat] = v1_new_4d_Observer0(z,u,y_m,par,ep)

% Observer based on the 4-state pipeline-riser model
% By: Esmaeil Jahanshahi
% May 2012, NTNU, Norway


% z1: Mass of gas in the pipeline (m_G1)
% z2: Mass of liquid in the pipeline (m_L1)
% z3: Pressure at top of riser (P_rt)
% z4: Mass of liquid in the riser (m_L2)

x1 = z(1); 
x2 = z(2); 
P2_t = z(3); 
x4 = z(4); 

u1 = u(1);
w_G_in = u(2);
w_L_in = u(3);

a = par.R*par.T2*par.rho_L/par.M_G;
b = par.rho_L*par.V2;


% rho_G1_norm = par.P1_norm*par.M_G/(par.R*par.T1);
% Alpha_L1_av = w_L_in*rho_G1_norm/(w_L_in*rho_G1_norm + w_G_in*par.rho_L);
h1ss = par.k_h*par.Alpha_L1_av*par.hc;
x2ss = par.V1*par.rho_L*par.Alpha_L1_av;
h1 = h1ss + sin(par.theta)*(x2 - x2ss)/(par.A1*(1-par.Alpha_L1_av)*par.rho_L);
h1 = max(h1,0);

V_G2 = par.V2 - x4/par.rho_L;


P1 = x1*par.R*par.T1/(par.M_G*(par.V1 - x2/par.rho_L));
rho_G1 = x1/(par.V1 - x2/par.rho_L);

Uslin = w_L_in/(par.A1*par.rho_L);
Re1=par.rho_L*Uslin*(2*par.r1)/par.visl;
Lambda1=0.0056 + 0.5*Re1^(-0.32);
Fric_pipe=0.5*par.Alpha_L1_av*Lambda1*par.rho_L*Uslin^2*par.L1/(2*par.r1);


%P2_t = x3*par.R*par.T2/(par.M_G*V_G2);
%P2_t = max(P2_t,par.P0);
x3 = P2_t*(b - x4)/a;

rho_G2 = x3/V_G2;
Alpha_L2_av = x4/(par.rho_L*par.V2);
rho_mix_av = (x3+x4)/par.V2;


Usl2 = w_L_in/(par.rho_L*par.A2);
Usg2 = w_G_in/(rho_G2*par.A2);
Um = Usl2+Usg2;

Re2=rho_mix_av*Um*(2*par.r2)/par.visl;
Lambda2=0.0056 + 0.5*Re2^(-0.32);
Fric_riser=0.5*Alpha_L2_av*Lambda2*rho_mix_av*Um^2*(par.L2+par.L3)/(2*par.r2);


A_g = (h1<par.hc)*(par.A1/par.hc^2)*(par.hc - h1)^2;

A_l = par.A1 - A_g;

P2_b = P2_t + rho_mix_av*par.g*par.L2 + Fric_riser;

w_G1 = par.K_g*A_g*sqrt(rho_G1*max(0,P1-Fric_pipe-P2_b));
w_L1 = par.K_o*A_l*sqrt(par.rho_L*max(0,P1-Fric_pipe+par.rho_L*par.g*h1-P2_b));

Alpha_Lb = 1 -  A_g/par.A1;

if(Alpha_Lb<=Alpha_L2_av)
    Alpha_Lb=Alpha_L2_av;
end

Alpha_Lt =  2*Alpha_L2_av - Alpha_Lb;


if(Alpha_Lt>Alpha_L2_av)
    Alpha_Lt = Alpha_L2_av;
elseif(Alpha_Lt<0)
    Alpha_Lt = 0;
end


Alpha_Lmt = Alpha_Lt*par.rho_L/(Alpha_Lt*par.rho_L + (1-Alpha_Lt)*rho_G2);
rho_t = Alpha_Lt*par.rho_L + (1-Alpha_Lt)*rho_G2; 

% if(z>0.95)
%     par.Cd = 0.95;
% end

%ORF = abs(1/(z^2*par.Cd^2) -1);
%w_mix_out = par.A2*sqrt(2*rho_t*max(0,P2_t-par.P0)/ORF);
w_mix_out = par.K_pc*u1*sqrt(rho_t*max(0,P2_t-par.P0));
%Q_out = w_mix_out/rho_t;

w_L_out = Alpha_Lmt*w_mix_out;
w_G_out = (1 - Alpha_Lmt)*w_mix_out;


dx1 = w_G_in - w_G1;
dx2 = w_L_in - w_L1;
dx3 = w_G1 - w_G_out;
dx4 = w_L1 - w_L_out;   

y_hat = P2_t/1e5; % Unit conversion from Pa to Bar

dz1 = dx1;
dz2 = dx2;
dz4 = dx4;
dz3 = (a*(b-x4)*dx3 + a*x3*dz4)/((b-x4)^2)+ (1/ep)*(y_m-y_hat);

zdot = [dz1;dz2;dz3;dz4];

xhat = [x1; x2; x3; x4];  
end