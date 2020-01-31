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
wG_in=u(2);
wL_in=u(3);

P1_ss = 70e5;
P2_start = 51e5;
VL2 = par.A2*par.L2/2;
x4_ss = VL2*par.rho_L;
VG2 = par.V2 - VL2;
rho_G2 = P2_start*par.M_G/(par.R*par.T2);
x3_ss = VG2*rho_G2;

%h1_ss = double(ss.h1(1));
rho_G1_n = P1_ss*par.M_G/(par.R*par.T1);
a_new = wL_in*rho_G1_n/(wL_in*rho_G1_n + wG_in*par.rho_L);
a_old = a_new - 0.01; % to make sure loop runs at least once

h1_old = 0.07;
h1_new = par.k_h*a_new*par.hc;

opt = optimset('Display','off','TolFun',1e-30,'TolX',1e-30,'MaxIter',1e4,'MaxFunEvals',1e10);
y0(1) = P1_ss/1e5;
while(abs(a_new-a_old)>1e-13)
    P1_ss = 1e5*y0(1);
    rho_G1 = P1_ss*par.M_G/(par.R*par.T1);
    par.Alpha_L1_av = wL_in*rho_G1/(wL_in*rho_G1 + wG_in*par.rho_L);
    x2_ss = par.rho_L*par.V1*a_new + (h1_old - h1_new)*par.A1*par.rho_L*(1-a_new)/sin(par.theta);
    x1_ss = (par.V1 - x2_ss/par.rho_L)*rho_G1;
    x0 = [x1_ss;x2_ss;x3_ss;x4_ss];
    [x0,fval,exitflag,output,jacobian] = fsolve(@v4_new_4d_model0,x0,opt,u,par);
    y0 = v4_new_4d_model(0,x0,u,'measurements',par);
    x3_ss = x0(3);
    x4_ss = x0(4);
    rho_G1 = 1e5*y0(1)*par.M_G/(par.R*par.T1);
    a_old = a_new;
    a_new = wL_in*rho_G1/(wL_in*rho_G1 + wG_in*par.rho_L);
    h1_old = h1_new;
    h1_new = par.k_h*a_new*par.hc;
end

str1 = ['Pressure in pipeline: ' num2str(y0(1)) ' [bar]' ];
disp(str1)
str2 = ['Pressure at top of riser: ' num2str(y0(2)) ' [bar]' ];
disp(str2)
str3 = ['Volumetric flow of choke: ' num2str(y0(3)) ' [L/s]' ];
disp(str3)
str4 = ['Mass flow of choke: ' num2str(y0(4)) ' [kg/s]' ];
disp(str4)

end