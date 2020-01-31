% *************** Finding Initial States *******************************
function [x0,y0,par] = v4_new_4d_initialize(u,par)

% clc
% clear all
% par = v1_new_4d_parameters();
% z = 0.6;
% wG_in = 0.36;
% wL_in = 8.64;

%  P1: Pressure in pipeline section
%  P2: Pressure at top of the riser
%  h1: Level of liquid in pipeline

%  x1: Mass of gas in the pipelien (m_G1)
%  x2: Mass of liquid in the pipeline (m_L1)
%  x3: Mass of gas in the riser (m_G2)
%  x4: Mass of liquid in the riser (m_L2)
z=u(1);
wG_in=u(2);
wL_in=u(3);
syms P1 x3 x4

h1 = 0.07;


%rho_G1 = P1*par.M_G/(par.R*par.T1);

P2_t = x3*par.R*par.T2/(par.M_G*(par.V2 - x4/par.rho_L));
%rho_G2 = x3/(par.V2 - x4/par.rho_L);
rho_G2 = P2_t*par.M_G/(par.R*par.T2);

Fric_pipe = 1.0164e+005;
Fric_riser = 3.6497e+004;

Alpha_L2_av = x4/(par.rho_L*par.V2);

rho_mix_av = (x3+x4)/par.V2;
P2_b = P2_t + rho_mix_av*par.g*par.L2 + Fric_riser;
A_g = par.A1 + par.A1*h1^2/par.hc^2 - 2*h1*par.A1/par.hc;
A_l = par.A1-A_g ;

%eq1 = wG_in^2 - (P1-Fric_pipe-P2_b)*rho_G1*(A_g*par.K_g)^2

eq2 = wL_in^2 - par.rho_L*(P1-Fric_pipe+par.rho_L*par.g*h1-P2_b)*(A_l*par.K_o)^2;


Alpha_Lt =  2*Alpha_L2_av - h1/par.hc;
rho_t = Alpha_Lt*par.rho_L + (1-Alpha_Lt)*rho_G2;
%alpha_Lm = H_LT*par.rho_L/(H_LT*par.rho_L + (1-H_LT)*rho_G2_T);

%if(z>0.9)
%    par.Cd = 0.9;
%end

%ORF = abs(1/(z^2*par.Cd^2) -1);

%eq3 = (wL_in+wG_in)^2 - (par.A2)^2*2*rho_t*(P2_t-par.P0)/ORF;
eq3 = (wL_in+wG_in)^2 - rho_t*(P2_t-par.P0)*(par.K_pc*z)^2;
%eq4 = simplify(wG_in^2 - rho_T*(P2_T-par.P_0)*((1-alpha_Lm)*par.K1*z)^2);
eq4 = wL_in/(wG_in + wL_in) - x4/(x4+x3);
%eq4 = rho_G2*wL_in/(rho_G2*wL_in + par.rho_L*wG_in) - Alpha_Lt;
%%
ss = solve(eq2,eq3,eq4);
%%
if(strcmp('7.4.0.287 (R2007a)',version))
    ind = 3;
else
    ind = 2;
end
P1_ss = real(double(ss.P1(ind)));
x3_ss = real(double(ss.x3(ind)));
x4_ss = real(double(ss.x4(ind)));
%h1_ss = double(ss.h1(1));
rho_G1_n = P1_ss*par.M_G/(par.R*par.T1);
a_new = wL_in*rho_G1_n/(wL_in*rho_G1_n + wG_in*par.rho_L);
a_old = a_new - 0.01; % to make sure loop runs at least once

h1_old = h1;
h1_new = par.k_h*a_new*par.hc;

opt = optimset('Display','off','TolFun',1e-20,'TolX',1e-20,'MaxIter',1e4,'MaxFunEvals',1e10);
y0(1) = P1_ss/1e5;
while(abs(a_new-a_old)>1e-13)
    par.P1_norm = 1e5*y0(1);
    rho_G1 = 1e5*y0(1)*par.M_G/(par.R*par.T1);
    par.Alpha_L1_av = wL_in*rho_G1/(wL_in*rho_G1 + wG_in*par.rho_L);
    x2_ss = par.rho_L*par.V1*a_new + (h1_old - h1_new)*par.A1*par.rho_L*(1-a_new)/sin(par.theta);
    x1_ss = (par.V1 - x2_ss/par.rho_L)*P1_ss*par.M_G/(par.R*par.T1);
    x0 = [x1_ss;x2_ss;x3_ss;x4_ss];
    [x0,fval,exitflag,output,jacobian] = fsolve(@v4_new_4d_model0,x0,opt,u,par);
    y0 = v4_new_4d_model(0,x0,u,'measurements',par);
    rho_G1 = 1e5*y0(1)*par.M_G/(par.R*par.T1);
    a_old = a_new;
    a_new = wL_in*rho_G1/(wL_in*rho_G1 + wG_in*par.rho_L);
    h1_old = h1_new;
    h1_new = par.k_h*a_new*par.hc;
end

str1 = ['Pressure in pipeline: ' num2str(y0(1)) ' bar' ];
disp(str1)
str2 = ['Pressure at top of riser: ' num2str(y0(2)) ' bar' ];
disp(str2)
str2 = ['Flow rate of choke: ' num2str(y0(3)) ' kg/s' ];
disp(str2)
end