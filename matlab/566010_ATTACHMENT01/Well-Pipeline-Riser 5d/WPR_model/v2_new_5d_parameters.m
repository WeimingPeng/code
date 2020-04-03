function par = v2_new_5d_parameters()
% ********** Constants and parameters for New simple model ***********
par.g = 9.81;                       % Gravity (m/s2)
par.R = 8314;                       % Gas constant (J/(K.Kmol))
par.rho_L = 832.2;                  % Liquid density (Kg/m3)
%par.rho_L = 832.2;
par.theta= pi/180;              % Feed pipe inclination (Rad) = 1 Deg
par.r_w = 0.06;                      % Raduis of well tubing (m)
par.D_w = 0.12;                      % Diameter of well tubing (m)
par.r_p = 0.06;                      % Raduis of pipe (m)
par.D_p = 0.12;                      % Diameter of pipe (m)
par.r_r = 0.05;                      % Raduis of riser (m)
par.D_r = 0.1;                       % Diameter of riser (m)
par.hc = 2*par.r_w/cos(par.theta);   % Critical liquid level (m)
par.Aw = pi*par.r_w^2;               % Cross section area of well tubing (m2)
par.Ap = pi*par.r_p^2;               % Cross section area of pipeline (m2)
par.Ar = pi*par.r_r^2;               % Cross section area of riser (m2)
par.Lw = 3000;                      % Height of oil well (m)
par.Lp = 4300;                      % Length of upstream pipe (m)
par.Lr = 300;                       % Height of riser
par.Lh = 100;                       % length of horozontal top section (m)
par.Tw = 337;                       % Well temprature (K)
par.Tp = 337;                       % Pipeline temprature (K)
%par.T1 = 300;
par.Tr = 298.3;                     % Riser temprature (K)
%par.T2=300;
par.M_Gw = 23.8;                    % Molecular weight of Gas (kg/kmol)
par.M_Gp = 22;
par.M_Gr = 23;

par.Pr = 320e5;                    % Reservior Pressure
par.P_s = 50.1e5;                    % Pressure after choke valve (Pa)

par.Vw = par.Aw*par.Lw;             % Total volume in upstream (m3)
par.Vp = par.Ap*par.Lp;    % Total volume in riser (m3)
par.Vr = par.Ar*(par.Lr+par.Lh);             % Total volume of well (m3)
par.PI = 2.75e-6;                     % Productivity Index of the well
% *************************************************************************
par.Pnormal = 101325;              % atmospheric conditions for Usg
Tnormal = 293.15;                  % Temprature in Normal conditions
RTnormal =par.R*Tnormal/par.M_Gr;   % copressibility in Normal Conditions
par.Rog_n = par.Pnormal/RTnormal;  % Gas density in normal conditions
par.visl= 1.4260e-4;
par.eps = 2.8e-5;                  % Roughness of pipe
par.Cd = 0.84;                        % Coef. of discharge
% *************** Finding Orifice Coefficients ****************************
% ******** Steady State Values Baesd on OLGA *************
u1 = 1;
u2 = 0.05;
P_bh = 284.902e5;
%w_r = par.PI*max(par.Pr - P_bh,0);
%w_Gr = 0;
P_wh = 74.7727e5;
%Alpha_Gwhm = 0.04175;
%w_Gwh = Alpha_Gwhm*w_r;
% -------------------------
w_in = 9.65183;
par.Alpha_Ginm = 0.0404106;
w_Gin = w_in*par.Alpha_Ginm;
w_Lin = w_in - w_Gin;
P_p = 74.6769e5;
%Alpha_L2_av = 0.717;
P_rt = 55.7673e5;
% --------------------------
w_Gout = w_Gin;
w_Lout = w_Lin;
w_mix = w_Gout + w_Lout;
% --------------------------
m_Gr = 58.7472;
m_Lr = 1693.99;
m_Gw = 379.205;
m_Lw = 23594.6;
m_Gp = 999.511;
m_Lp = 26511.6;

% for Z2=0.1
% x1 = 23731
% x2 = 897.54
% x3 = 26935.6
% x3 = 59.274
% x4 = 1594.6
% ********************************************************
par.w_nom = 13;

V_Gp = par.Vp - m_Lp/par.rho_L;
rho_Gp = m_Gp/V_Gp;
%P_p = par.R*par.Tp*rho_Gp/par.M_Gp; % Pressure at inlet of pipeline
%par.M_Gp = par.R*par.Tp*rho_Gp/P_p;
Alpha_Gmp = 0.035; %par.Alpha_Ginm;
par.Alpha_Lp_av = rho_Gp*(1-Alpha_Gmp)/(rho_Gp*(1-Alpha_Gmp)+par.rho_L*Alpha_Gmp);

Alpha_Gmw_av = m_Gw/(m_Gw+m_Lw);

Alpha_Lw_av = m_Lw/(par.rho_L*par.Vw);
par.k_a = 2*Alpha_Gmw_av/Alpha_Gmp;

V_Gw = par.Vw - m_Lw/par.rho_L;
rho_Gw = m_Gw/V_Gw;  % Density of gas at well-head

%par.M_Gw = rho_Gw*par.R*par.Tp/P_wh;

rho_mix_w = (m_Gw+m_Lw)/par.Vw;

Alpha_Gwmt = Alpha_Gmp;
Alpha_Lwt = (1-Alpha_Gwmt)*rho_Gw/((1-Alpha_Gwmt)*rho_Gw + Alpha_Gwmt*par.rho_L);
rho_mixwt = Alpha_Lwt*rho_Gw + (1-Alpha_Lwt)*par.rho_L;


P_bh = par.Pr - w_in/par.PI;

Fric = P_bh - P_wh - rho_mix_w*par.g*par.Lw;
Uslw = par.w_nom/(par.Aw*par.rho_L);

Re_w=par.rho_L*Uslw*par.D_w/par.visl;

% temp = -1.8*log10((par.eps/par.D_w/3.7)^1.11+6.9/Re_w);
% Lambda_w=(1/temp)^2;

Lambda_w=0.0056 + 0.5*Re_w^(-0.32);
Fric_w=0.5*Alpha_Lw_av*Lambda_w*par.rho_L*Uslw^2*par.Lw/(2*par.r_w);

P_bh = P_wh + rho_mix_w*par.g*par.Lw + Fric_w;


% rho_G1_norm = par.P1_norm*par.M_Gp/(par.R*par.T1);
% Alpha_L1_av = w_L_in*rho_G1_norm/(w_L_in*rho_G1_norm + w_G_in*par.rho_L);
par.k_h = 0.6;

hss = par.k_h*par.Alpha_Lp_av*par.hc;
h = hss;


% if(par.Cd==1 && u1>0.95)
%     u1 = 0.95;
% end
%ORF = abs(1/(u1^2*par.Cd^2) -1);
%w_mix_wh = par.Ap*sqrt(2*rho_mixwt*max(0,P_wh-P_p)/ORF); 
    
par.K1 = 2.5*w_mix/(u1*sqrt(rho_mixwt*max(0,P_wh-P_p))); % Flow rate form well to pipeline
%Alpha_Gmp = rho_Gp*(1-par.Alpha_Lp_av)/(rho_Gp*(1-par.Alpha_Lp_av)+par.rho_L*par.Alpha_Lp_av);

Uslp = w_Lin/(par.Ap*par.rho_L);
Re_p=par.rho_L*Uslp*(2*par.r_p)/par.visl;
%Lambda_p=0.0056 + 0.5*Re_p^(-0.32);

temp = -1.8*log10((par.eps/par.D_p/3.7)^1.11+6.9/Re_p);
Lambda_p=(1/temp)^2;

Fric_p=0.5*par.Alpha_Lp_av*Lambda_p*par.rho_L*Uslp^2*par.Lp/(2*par.r_p);

V_Gr = par.Vr - m_Lr/par.rho_L;

% P_rt = max(P_rt,par.P_s);
% m_Gr = P_rt*(b3 - m_Lr)/a3;

%par.M_Gr = m_Gr*par.R*par.Tr/(P_rt*V_Gr);

rho_Gr = m_Gr/V_Gr;
Alpha_Lr_av = m_Lr/(par.rho_L*par.Vr);

A_g = (h<par.hc)*(par.Ap/par.hc^2)*(par.hc - h)^2;

A_l = par.Ap - A_g;

Alpha_Lrb = 1 -  A_g/par.Ap;

if(Alpha_Lrb<=Alpha_Lr_av)
    Alpha_Lrb=Alpha_Lr_av;
end

Alpha_Lrt =  2*Alpha_Lr_av - Alpha_Lrb;


if(Alpha_Lrt>Alpha_Lr_av)
    Alpha_Lrt = Alpha_Lr_av;
elseif(Alpha_Lrt<0)
    Alpha_Lrt = 0;
end


%rho_mixr_av1 = (m_Gr+m_Lr-Alpha_Lrt*par.Lh*par.Ar*par.rho_L-(1-Alpha_Lrt)*par.Lh*par.Ar*rho_Gr)/(par.Ar*par.Lr)
rho_mixr_av = (m_Gr+m_Lr)/par.Vr;

Uslr = w_Lin/(par.rho_L*par.Ar);
Usgr = w_Gin/(rho_Gr*par.Ar);
Umr = Uslr+Usgr;

Re_r=rho_mixr_av*Umr*(2*par.r_r)/par.visl;
%Lambda_r=0.0056 + 0.5*Re_r^(-0.32);

temp = -1.8*log10((par.eps/par.D_r/3.7)^1.11+6.9/Re_r);
Lambda_r=(1/temp)^2;

Fric_r=0.5*Alpha_Lr_av*Lambda_r*Umr^2*(par.Lr+par.Lh)/(2*par.r_r);


DP_r = rho_mixr_av*par.g*par.Lr + Fric_r;
P_rb = P_rt + DP_r;


par.K_g = 1.2*w_Gin/(A_g*sqrt(rho_Gp*max(0,P_p-P_rb)));
par.K_o = 1.7*w_Lin/(A_l*sqrt(par.rho_L*abs(P_p-Fric_p+par.rho_L*par.g*h-P_rb)));

rho_mixrt = Alpha_Lrt*par.rho_L + (1-Alpha_Lrt)*rho_Gr; 

% if(z>0.95)
%     par.Cd = 0.95;
% end

%ORF = abs(1/(z^2*par.Cd^2) -1);
%w_mix_out = par.A2*sqrt(2*rho_t*max(0,P2_t-par.P0)/ORF);
par.K2 = 1*w_mix/(u2*sqrt(rho_mixrt*max(0,P_rt-par.P_s)));


end