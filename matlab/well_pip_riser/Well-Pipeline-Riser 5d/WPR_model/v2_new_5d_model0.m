function [xdot,y] = v2_new_5d_model0(x,u,par)

% Simplified model for riser slugging
% Oil well added as inlet boundary condition
% By: Esmaeil Jahanshahi
% August 2011, NTNU, Norway

% x1: Total mass in the well (m_w)
% x2: Mass of gas in the pipeline (m_Gp)
% x3: Mass of liquid in the pipeline (m_Lp)
% x4: Mass of gas in the riser (m_Gr)
% x5: Mass of liquid in the riser (m_Lr)

m_w = x(1); 
m_Gp = x(2); 
m_Lp = x(3); 
m_Gr = x(4); 
m_Lr = x(5); 

u1 = u(1);
u2 = u(2);
d1 = u(3);
d2 = u(4);

V_Gp = par.Vp - m_Lp/par.rho_L;
rho_Gp = m_Gp/V_Gp;
P_p = par.R*par.Tp*rho_Gp/par.M_Gp; % Pressure at inlet of pipeline

Alpha_Gmp = rho_Gp*(1-par.Alpha_Lp_av)/(rho_Gp*(1-par.Alpha_Lp_av)+par.rho_L*par.Alpha_Lp_av);
Alpha_Gwmt = Alpha_Gmp;

Alpha_Gmw_av = par.k_a*Alpha_Gwmt/2;

m_Gw = Alpha_Gmw_av*m_w;
m_Lw = (1-Alpha_Gmw_av)*m_w;

Alpha_Lw_av = m_Lw/(par.rho_L*par.Vw);

V_Gw = par.Vw - m_Lw/par.rho_L;

P_wh = m_Gw*par.R*par.Tp/(par.M_Gw*V_Gw); % Pressure at well-head
rho_Gw = m_Gw/V_Gw;  % Density of gas at well-head
rho_mix_w = m_w/par.Vw;

Uslw = par.w_nom/(par.Aw*par.rho_L);

Re_w=par.rho_L*Uslw*par.D_w/par.visl;

temp = -1.8*log10((par.eps/par.D_w/3.7)^1.11+6.9/Re_w);
Lambda_w=(1/temp)^2;

%Lambda_w=0.0056 + 0.5*Re_w^(-0.32);
Fric_w=0.5*Alpha_Lw_av*Lambda_w*par.rho_L*Uslw^2*par.Lw/(2*par.r_w);
P_bh = P_wh + rho_mix_w*par.g*par.Lw + Fric_w;
w_r = par.PI*max(par.Pr - P_bh,0);
%Alpha_Gwmt = par.b - par.a*P_wh;
%Alpha_Gwmt = 4*m_Gw/(m_Gw+m_Lw);
%Alpha_Gwmt = 0.04;

% if(Alpha_Lwt>Alpha_Lw_av)
%     Alpha_Lwt = Alpha_Lw_av;
% elseif(Alpha_Lwt<0)
%     Alpha_Lwt = 0;
% end

Alpha_Lwt = (1-Alpha_Gwmt)*rho_Gw/(Alpha_Gwmt*par.rho_L + (1-Alpha_Gwmt)*rho_Gw);

rho_mixwt = (1-Alpha_Lwt)*rho_Gw + Alpha_Lwt*par.rho_L;


% rho_G1_norm = par.P1_norm*par.M_Gp/(par.R*par.T1);
% Alpha_L1_av = w_L_in*rho_G1_norm/(w_L_in*rho_G1_norm + w_G_in*par.rho_L);


% if(par.Cd==1 && u1>0.95)
%     u1 = 0.95;
% end
%ORF = abs(1/(u1^2*par.Cd^2) -1);
%w_mix_wh = par.Ap*sqrt(2*rho_mixwt*max(0,P_wh-P_p)/ORF); 
    
w_mix_wh = par.K1*u1*sqrt(rho_mixwt*max(0,P_wh-P_p)); % Flow rate form well to pipeline

w_Gwh = w_mix_wh*Alpha_Gmp;
w_Gin = w_Gwh + d1;

w_Lwh = w_mix_wh*(1 - Alpha_Gmp);
w_Lin = w_Lwh + d2;

hss = par.k_h*par.Alpha_Lp_av*par.hc;
m_Lpss = par.Vp*par.rho_L*par.Alpha_Lp_av;
K = sin(par.theta)/(par.Ap*(1-par.Alpha_Lp_av)*par.rho_L);
h = max(hss + K*(m_Lp - m_Lpss),0);

Uslp = w_Lin/(par.Ap*par.rho_L);
Re_p=par.rho_L*Uslp*(2*par.r_p)/par.visl;
%Lambda_p=0.0056 + 0.5*Re_p^(-0.32);

temp = -1.8*log10((par.eps/par.D_p/3.7)^1.11+6.9/Re_p);
Lambda_p=(1/temp)^2;

Fric_p=0.5*par.Alpha_Lp_av*Lambda_p*par.rho_L*Uslp^2*par.Lp/(2*par.r_p);

V_Gr = par.Vr - m_Lr/par.rho_L;

% P_rt = max(P_rt,par.P_s);
% m_Gr = P_rt*(b3 - m_Lr)/a3;

P_rt = m_Gr*par.R*par.Tr/(par.M_Gr*V_Gr);
P_rt = max(P_rt,par.P_s);

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


w_G_lp = par.K_g*A_g*sqrt(rho_Gp*max(0,P_p-P_rb));
w_L_lp = (P_p-Fric_p+par.rho_L*par.g*h-P_rb>0)*...
    par.K_o*A_l*sqrt(par.rho_L*abs(P_p-Fric_p+par.rho_L*par.g*h-P_rb));


Alpha_Lrmt = Alpha_Lrt*par.rho_L/(Alpha_Lrt*par.rho_L + (1-Alpha_Lrt)*rho_Gr);
rho_mixrt = Alpha_Lrt*par.rho_L + (1-Alpha_Lrt)*rho_Gr; 

% if(z>0.95)
%     par.Cd = 0.95;
% end

%ORF = abs(1/(z^2*par.Cd^2) -1);
%w_mix_out = par.A2*sqrt(2*rho_t*max(0,P2_t-par.P0)/ORF);
w_mix_out = par.K2*u2*sqrt(rho_mixrt*max(0,P_rt-par.P_s));
Q_out = 1000*w_mix_out/rho_mixrt;
w_L_out = Alpha_Lrmt*w_mix_out;
w_G_out = (1 - Alpha_Lrmt)*w_mix_out;

dx1 = w_r - w_mix_wh;
dx2 = w_Gin - w_G_lp;
dx3 = w_Lin - w_L_lp;
dx4 = w_G_lp - w_G_out;
dx5 = w_L_lp - w_L_out;


xdot = [dx1;dx2;dx3;dx4;dx5];

y = [P_bh     %Y1 [pa]
     P_wh     %Y2 [pa]
     w_mix_wh %Y3 [kg/s]
     P_p      %Y4 [Pa]
     P_rb     %Y5 [pa]
     DP_r     %Y6 [Pa]
     P_rt     %Y7 [kg/m3]
     Q_out    %Y8 [Litre/s]
     w_mix_out %Y9 [kg/s]
     rho_mixrt%Y10[kg/m3]
     Alpha_Lrt%Y11[-]
     w_Gin    
     w_Lin];

end