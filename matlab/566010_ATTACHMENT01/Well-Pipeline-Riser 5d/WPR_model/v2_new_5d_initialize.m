% *************** Finding Initial States *******************************
function [x0,y0,par] = v2_new_5d_initialize(u,par)

% clc
% clear all
% par = v1_new_6d_parameters();
% z1 = 1;
% z2 = 0.04;
% d1 = 0;
% d2 = 0;
% u = [z1;z2;d1;d2];


% x1: Mass of gas in the well (m_Gw)
% x2: Mass of gas in the pipelien (m_G1)
% x3: Mass of liquid in the pipeline (m_L1)
% x4: Mass of gas in the riser (m_G2)
% x5: Mass of liquid in the riser (m_L2)

% x1 = 403.148;
% x2 = 883.54;
% x3 = 26967.3;
% x4 = 59.2237;
% x5 = 1577.78;
if (u(2)<0.1)
    P_wh = 79.2653e5;
    P_bh = 291.292e5;
    P_p_ss = 79.1709e5;
    P_rt_ss = 60.6255e5;
elseif(0.06<u(2)<0.6)
    P_wh = 7019947;
    P_bh = 27854510;
    P_p_ss = 7011254;
    P_rt_ss = 5069538;
else
    P_wh =6980175;
    P_bh = 27799650;
    P_p_ss = 6971484;
    P_rt_ss = 5018719;
end

rho_Gp = P_p_ss*par.M_Gp/(par.R*par.Tp);
Alpha_Gmp = rho_Gp*(1-par.Alpha_Lp_av)/(rho_Gp*(1-par.Alpha_Lp_av)+par.rho_L*par.Alpha_Lp_av);
Alpha_Gmw_av = par.k_a*Alpha_Gmp/2;

a = Alpha_Gmw_av/(1-Alpha_Gmw_av);

m_Lw_ss = P_p_ss*par.M_Gw*par.Vw/(P_p_ss*par.M_Gw/par.rho_L + a*par.R*par.Tp); 
m_Gw_ss = a*m_Lw_ss;

   

VLr = par.Ar*par.Lr/2;
m_Lr_ss = VLr*par.rho_L;
VGr = par.Vr - VLr;
rho_Gr = P_rt_ss*par.M_Gr/(par.R*par.Tr);
m_Gr_ss = VGr*rho_Gr;

m_Lp_ss = par.rho_L*par.Vp*par.Alpha_Lp_av + (0.01)*par.Ap*par.rho_L*(1-par.Alpha_Lp_av)/sin(par.theta);
rho_Gp = P_p_ss*par.M_Gp/(par.R*par.Tp);
m_Gp_ss = (par.Vp - m_Lp_ss/par.rho_L)*rho_Gp;

x0 = [m_Gw_ss+m_Lw_ss m_Gp_ss m_Lp_ss m_Gr_ss m_Lr_ss]';

opt = optimset('Display','off','TolFun',1e-50,'TolX',1e-50,'MaxIter',1e4,'MaxFunEvals',1e10);

[x0,fval,exitflag,output,jacobian] = fsolve(@v2_new_5d_model0,x0,opt,u,par);
y0 = v2_new_5d_model(0,x0,u,'measurements',par);

% opt = optimset('Display','off','TolFun',1e-50,'TolX',1e-50,'MaxIter',1e4,'MaxFunEvals',1e10);
%
% x0 = [x1;x2;x3;x4;x5];
% [x0,fval,exitflag,output,jacobian] = fsolve(@v1_new_5d_model0,x0,opt,u,par)
% y0 = v1_new_5d_model(0,x0,u,'measurements',par);
%par.Alpha_Lp_av = a_new;
str1 = ['Pressure at bottom-hole (P_bh):' num2str(y0(1)) ' [bar]' ];
disp(str1)
str2 = ['Pressure at well-head (P_wh):' num2str(y0(2)) ' [bar]' ];
disp(str2)
str3 = ['Inlet mass flow (w_in):' num2str(y0(3)) ' [kg/s]' ];
disp(str3)
str4 = ['Pressure at inlet of pipeline (P_in):' num2str(y0(4)) ' [bar]' ];
disp(str4)
str5 = ['Pressure at top of riser (P2_t):' num2str(y0(7)) ' [bar]' ];
disp(str5)
str6 = ['Mass flow of choke (w_out):' num2str(y0(9)) ' [kg/s]' ];
disp(str6)
end