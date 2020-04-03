% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Linerized model of final version New model. Function version.
%
% x1=mg1
% x2=ml1
% x3=mg2
% x4=ml2
% u1=z
% u2=w_G_in
% u3=w_L_in
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Written by:  Knute Meland
% Date:         27.10.2010
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function [A,B,C,D]= v4_new_4d_linmod_knut(x0,y0,u,par)

%%
u1 = u(1);
u2 = u(2);
u3 = u(3);

x1 = x0(1);
x2 = x0(2);
x3 = x0(3);
x4 = x0(4);

P1_nom=y0(1)*10^5;

%%
g=par.g;
R=par.R;
rho_L=par.rho_L;
theta=par.theta; 
r1=par.r1;
r2=par.r2; 
hc=par.hc;
A1=par.A1; 
A2=par.A2;
L1=par.L1; 
L2=par.L2;
L3=par.L3;
T1=par.T1;
T2=par.T2;
M_G=par.M_G;
P0=par.P0;
V2=par.V2; 
V1=par.V1;

K_pc =par.K_pc;
K_g =par.K_g;
K_L=par.K_o; 

K_h=par.k_h; %correction factor to finetune model
my=par.visl; %friction factor 



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating differensials:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rho_G1_av=P1_nom*M_G/(R*T1);

%%
a_L1_av=(rho_G1_av*u3)/(rho_G1_av*u3 + rho_L*u2);

da_L1_av__dx1 = 0;
da_L1_av__dx2 = 0;
da_L1_av__dx3 = 0;
da_L1_av__dx4 = 0;
da_L1_av__du1 = 0;
da_L1_av__du2 = 0; %-(rho_L*rho_G1_av*u3)/(rho_G1_av*u3 +rho_L*u2)^2;
da_L1_av__du3 = 0; %(rho_G1_av*rho_L*u2)/(rho_G1_av*u3 +rho_L*u2)^2;
%%
h1_av = K_h*hc*a_L1_av;

dh1_av__dx1 = 0;
dh1_av__dx2 = 0;
dh1_av__dx3 = 0;
dh1_av__dx4 = 0;
dh1_av__du1 = 0;
dh1_av__du2 = K_h*hc*da_L1_av__du2;
dh1_av__du3 = K_h*hc*da_L1_av__du3;
%%
h1 = h1_av + sin(theta)*(x2 - V1*rho_L*a_L1_av)/(A1*(1-a_L1_av)*rho_L);

dh1__dx1 = 0;
dh1__dx2 = sin(theta)/(A1*(1-a_L1_av)*par.rho_L);
dh1__dx3 = 0;
dh1__dx4 = 0;
dh1__du1 = 0;
dh1__du2 = dh1_av__du2 ...
    +sin(theta)*x2*da_L1_av__du2/(A1*rho_L*(1-a_L1_av)^2)...
    -(V1*sin(theta)/A1)*da_L1_av__du2/(1-a_L1_av)^2;
dh1__du3 = dh1_av__du3 ...
    +sin(theta)*x2*da_L1_av__du3/(A1*rho_L*(1-a_L1_av)^2)...
    -(V1*sin(theta)/A1)*da_L1_av__du3/(1-a_L1_av)^2;

%%
V_G1 = V1 -x2/rho_L;

dV_G1__dx1 = 0;
dV_G1__dx2 = -1/rho_L;
dV_G1__dx3 = 0;
dV_G1__dx4 = 0;
dV_G1__du1 = 0;
dV_G1__du2 = 0;
dV_G1__du3 = 0;

%%
rho_G1 = x1/V_G1;

drho_G1__dx1 = 1/V_G1;
drho_G1__dx2 = -x1*dV_G1__dx2/V_G1^2;
drho_G1__dx3 = 0;
drho_G1__dx4 = 0;
drho_G1__du1 = 0;
drho_G1__du2 = 0;
drho_G1__du3 = 0;

%%
P1 = rho_G1*R*T1/M_G;

dP1__dx1 = R*T1*drho_G1__dx1 / M_G;
dP1__dx2 = R*T1*drho_G1__dx2 / M_G;
dP1__dx3 = 0;
dP1__dx4 = 0;
dP1__du1 = 0;
dP1__du2 = 0;
dP1__du3 = 0;

%%
U_slin_av = u3/(pi*r1^2*rho_L);

dU_slin_av__dx1 = 0;
dU_slin_av__dx2 = 0;
dU_slin_av__dx3 = 0;
dU_slin_av__dx4 = 0;
dU_slin_av__du1 = 0;
dU_slin_av__du2 = 0;
dU_slin_av__du3 = 1 / (pi*r1^2*rho_L );

%%
Re_p = 2*rho_L*U_slin_av*r1 / my ;

dRe_p__dx1 = 0;
dRe_p__dx2 = 0;
dRe_p__dx3 = 0;
dRe_p__dx4 = 0;
dRe_p__du1 = 0;
dRe_p__du2 = 0;
dRe_p__du3 = 2*rho_L*r1* dU_slin_av__du3 / my;

%%
L_p = 0.0056 +0.5*Re_p^(-0.32);

dL_p__dx1 = 0;
dL_p__dx2 = 0;
dL_p__dx3 = 0;
dL_p__dx4 = 0;
dL_p__du1 = 0;
dL_p__du2 = 0;
dL_p__du3 = -0.5*0.32*Re_p^(-1.32)*dRe_p__du3;

%%
DP_fp = a_L1_av*L_p*rho_L*U_slin_av^2*L1/(4*r1);

dDP_fp__dx1 = 0;
dDP_fp__dx2 = 0;
dDP_fp__dx3 = 0;
dDP_fp__dx4 = 0;
dDP_fp__du1 = 0;
dDP_fp__du2 = L_p*rho_L*U_slin_av^2*L1*da_L1_av__du2 / (4*r1);
dDP_fp__du3 = ( L_p*rho_L*U_slin_av^2*L1*da_L1_av__du3  ...
    +a_L1_av*rho_L*U_slin_av^2*L1*dL_p__du3 ...
    +2*a_L1_av*rho_L*L_p*U_slin_av*L1*dU_slin_av__du3 ) / (4*r1);

%%
V_G2 = V2 - x4 / rho_L ; 

dV_G2__dx1 = 0;
dV_G2__dx2 = 0;
dV_G2__dx3 = 0;
dV_G2__dx4 = -1/rho_L ;
dV_G2__du1 = 0;
dV_G2__du2 = 0;
dV_G2__du3 = 0;

%%
rho_G2 = x3 / V_G2 ; 

drho_G2__dx1 = 0;
drho_G2__dx2 = 0;
drho_G2__dx3 = 1/V_G2;
drho_G2__dx4 = -x3*dV_G2__dx4 / V_G2^2;
drho_G2__du1 = 0;
drho_G2__du2 = 0;
drho_G2__du3 = 0;

%%
P2 = rho_G2*R*T2 / M_G; 

dP2__dx1 = 0;
dP2__dx2 = 0;
dP2__dx3 = R*T2*drho_G2__dx3/M_G;
dP2__dx4 = R*T2*drho_G2__dx4/M_G;
dP2__du1 = 0;
dP2__du2 = 0;
dP2__du3 = 0;

%%
Alpha_L2_av = x4/(rho_L*V2);

dAlpha_L2__dx1 = 0; 
dAlpha_L2__dx2 = 0;
dAlpha_L2__dx3 = 0;
dAlpha_L2__dx4 = 1/(rho_L*V2);
dAlpha_L2__du1 = 0;
dAlpha_L2__du2 = 0;
dAlpha_L2__du3 = 0;

%%
rho_m_av = (x3+x4)/V2;

drho_m_av__dx1 = 0; 
drho_m_av__dx2 = 0;
drho_m_av__dx3 = 1/V2;
drho_m_av__dx4 = 1/V2;
drho_m_av__du1 = 0; 
drho_m_av__du2 = 0; 
drho_m_av__du3 = 0; 

%%
U_sl2_av=u3/(rho_L*pi*r2^2);
U_sg2_av=u2 /(rho_G2*pi*r2^2);

%%
U_m_av = U_sl2_av +U_sg2_av;

dU_m_av__dx1 = 0;
dU_m_av__dx2 = 0;
dU_m_av__dx3 = -U_sg2_av * drho_G2__dx3 / rho_G2;
dU_m_av__dx4 = -U_sg2_av * drho_G2__dx4 / rho_G2;
dU_m_av__du1 = 0;
dU_m_av__du2 = U_sg2_av / u2 ; 
dU_m_av__du3 = U_sl2_av / u3 ; 

%%
Re_r= 2*rho_m_av *U_m_av* r2 / my;

dRe_r__dx1 = 0;
dRe_r__dx2 = 0;
dRe_r__dx3 = 2*rho_m_av * r2* dU_m_av__dx3 / my  ...
    + 2 *U_m_av* r2* drho_m_av__dx3 / my;
dRe_r__dx4 = 2*rho_m_av * r2* dU_m_av__dx4 / my  ...
    + 2 *U_m_av* r2* drho_m_av__dx4 / my;
dRe_r__du1 = 0;
dRe_r__du2 = 2*rho_m_av * r2* dU_m_av__du2 / my;
dRe_r__du3 = 2*rho_m_av * r2* dU_m_av__du3 / my;

%%
L_r = 0.0056 + 0.5*Re_r^(-0.32);

dL_r__dx1 = 0;
dL_r__dx2 = 0;
dL_r__dx3 = -0.5*0.32*Re_r^(-1.32)*dRe_r__dx3;
dL_r__dx4 = -0.5*0.32*Re_r^(-1.32)*dRe_r__dx4;
dL_r__du1 = 0;
dL_r__du2 = -0.5*0.32*Re_r^(-1.32)*dRe_r__du2;
dL_r__du3 = -0.5*0.32*Re_r^(-1.32)*dRe_r__du3;

%% 
DP_fr = Alpha_L2_av * L_r *rho_m_av *U_m_av^2 * (L2 + L3) / (4*r2);

dDP_fr__dx1 = 0;
dDP_fr__dx2 = 0;
dDP_fr__dx3 = ( 2*Alpha_L2_av * L_r *rho_m_av *U_m_av * (L2 + L3)*dU_m_av__dx3 ...
    + Alpha_L2_av * L_r *U_m_av^2 * (L2 + L3) * drho_m_av__dx3 ...
    +Alpha_L2_av *rho_m_av *U_m_av^2 * (L2 + L3) *dL_r__dx3 ) / (4*r2);
dDP_fr__dx4 = ( 2*Alpha_L2_av * L_r *rho_m_av *U_m_av * (L2 + L3)*dU_m_av__dx4 ...
    + Alpha_L2_av * L_r *U_m_av^2 * (L2 + L3) * drho_m_av__dx4 ...
    +Alpha_L2_av *rho_m_av *U_m_av^2 * (L2 + L3) *dL_r__dx4 ...
    +L_r *rho_m_av *U_m_av^2 * (L2 + L3)*dAlpha_L2__dx4 )/ (4*r2);
dDP_fr__du1 = 0;
dDP_fr__du2 = ( 2*Alpha_L2_av * L_r *rho_m_av *U_m_av * (L2 + L3)*dU_m_av__du2 ...
    +Alpha_L2_av *rho_m_av *U_m_av^2 * (L2 + L3) *dL_r__du2 ) / (4*r2);
dDP_fr__du3 = ( 2*Alpha_L2_av * L_r *rho_m_av *U_m_av * (L2 + L3)*dU_m_av__du3 ...
    +Alpha_L2_av *rho_m_av *U_m_av^2 * (L2 + L3) *dL_r__du3 ) / (4*r2);

%%
DP_G = P1 -DP_fp - P2 - rho_m_av * g * L2 - DP_fr;

dDP_G__dx1 = dP1__dx1;
dDP_G__dx2 = dP1__dx2; 
dDP_G__dx3 = - dP2__dx3 - g*L2*drho_m_av__dx3 - dDP_fr__dx3;
dDP_G__dx4 = - dP2__dx4 - g*L2*drho_m_av__dx4 - dDP_fr__dx4;
dDP_G__du1 = 0;
dDP_G__du2 = -dDP_fp__du2  - dDP_fr__du2;
dDP_G__du3 = -dDP_fp__du3  - dDP_fr__du3;

%%
DP_L = P1 -DP_fp +rho_L*g*h1 - P2 -rho_m_av*g *L2 - DP_fr;

dDP_L__dx1 = dP1__dx1;
dDP_L__dx2 = dP1__dx2 + rho_L* g* dh1__dx2;
dDP_L__dx3 = -dP2__dx3 - g* L2 * drho_m_av__dx3 - dDP_fr__dx3;
dDP_L__dx4 = -dP2__dx4 - g* L2 * drho_m_av__dx4 - dDP_fr__dx4;
dDP_L__du1 = 0;
dDP_L__du2 = -dDP_fp__du2 + rho_L*g*dh1__du2 - dDP_fr__du2;
dDP_L__du3 = -dDP_fp__du3 +rho_L*g*dh1__du3 - dDP_fr__du3;

%%
if h1<hc 
    A_G = (A1/hc^2)*(hc - h1)^2; 

    dA_G__dx1 = 0;
    dA_G__dx2 = (-2*A1/hc^2)*(hc - h1)*dh1__dx2;
    dA_G__dx3 = 0;
    dA_G__dx4 = 0;
    dA_G__du1 = 0;
    dA_G__du2 = (-2*A1/hc^2)*(hc - h1)*dh1__du2;
    dA_G__du3 = (-2*A1/hc^2)*(hc - h1)*dh1__du3;
    
else
    A_G=0;
    
    dA_G__dx1 = 0;
    dA_G__dx2 = 0;
    dA_G__dx3 = 0;
    dA_G__dx4 = 0;
    dA_G__du1 = 0;
    dA_G__du2 = 0;
    dA_G__du3 = 0;
    
end

%%
A_L = A1 - A_G;

dA_L__dx1 = 0;
dA_L__dx2 = -dA_G__dx2;
dA_L__dx3 = 0;
dA_L__dx4 = 0;
dA_L__du1 = 0;
dA_L__du2 = -dA_G__du2;
dA_L__du3 = -dA_G__du3;

%%
if h1<hc
    if DP_G<0
        disp('The flowdirection is wrong. Check parameters');
        return
    end
    w_G_lp = K_g*A_G*sqrt(rho_G1*DP_G); 

    dw_G_lp__dx1 = K_g*A_G*(DP_G*drho_G1__dx1+ rho_G1*dDP_G__dx1)/(2*sqrt(rho_G1*DP_G)); 
    dw_G_lp__dx2 = K_g*dA_G__dx2*sqrt(rho_G1*DP_G)...
                + K_g*A_G*(DP_G*drho_G1__dx2 +rho_G1*dDP_G__dx2)/(2*sqrt(rho_G1*DP_G));
    dw_G_lp__dx3 = K_g*A_G*rho_G1*dDP_G__dx3/(2*sqrt(rho_G1*DP_G));
    dw_G_lp__dx4 = K_g*A_G*rho_G1*dDP_G__dx4/(2*sqrt(rho_G1*DP_G));
    dw_G_lp__du1 = 0;
    dw_G_lp__du2 = K_g*A_G*rho_G1*dDP_G__du2/(2*sqrt(rho_G1*DP_G))...
        +K_g*dA_G__du2*sqrt(rho_G1*DP_G);
    dw_G_lp__du3 = K_g*A_G*rho_G1*dDP_G__du3/(2*sqrt(rho_G1*DP_G))...
        +K_g*dA_G__du3*sqrt(rho_G1*DP_G);
else
    w_G_lp = 0;

    dw_G_lp__dx1 = 0;
    dw_G_lp__dx2 = 0;
    dw_G_lp__dx3 = 0;
    dw_G_lp__dx4 = 0;
    dw_G_lp__du1 = 0;
    dw_G_lp__du2 = 0;
    dw_G_lp__du3 = 0;
end

%%
if DP_L<0
    disp('The flowdirection is wrong. Check parameters 2');
    return
end
    
w_L_lp = K_L*A_L*sqrt(rho_L*DP_L);

dw_L_lp__dx1 = K_L*A_L*rho_L*dDP_L__dx1/(2*sqrt(rho_L*DP_L)); 
dw_L_lp__dx2 = K_L*A_L*rho_L*dDP_L__dx2/(2*sqrt(rho_L*DP_L)) ...
    +K_L*dA_L__dx2*sqrt(rho_L*DP_L);
dw_L_lp__dx3 = K_L*A_L*rho_L*dDP_L__dx3/(2*sqrt(rho_L*DP_L));
dw_L_lp__dx4 = K_L*A_L*rho_L*dDP_L__dx4/(2*sqrt(rho_L*DP_L));
dw_L_lp__du1 = 0;
dw_L_lp__du2 = K_L*A_L*rho_L*dDP_L__du2/(2*sqrt(rho_L*DP_L)) ...
    +K_L*dA_L__du2*sqrt(rho_L*DP_L);
dw_L_lp__du3 = K_L*A_L*rho_L*dDP_L__du3/(2*sqrt(rho_L*DP_L)) ...
    +K_L*dA_L__du3*sqrt(rho_L*DP_L);

%%
Alpha_Lt = 2*x4 / (V2*rho_L) - A_L/A1 ;

dAlpha_Lt__dx1 = 0;
dAlpha_Lt__dx2 = -dA_L__dx2/A1;
dAlpha_Lt__dx3 = 0;
dAlpha_Lt__dx4 = 2/(V2*rho_L);
dAlpha_Lt__du1 = 0;
dAlpha_Lt__du2 = -dA_L__du2/A1;
dAlpha_Lt__du3 = -dA_L__du3/A1;

%%
rho_t=Alpha_Lt*rho_L +(1-Alpha_Lt)*rho_G2;

drho_t__dx1 = 0;
drho_t__dx2 = (rho_L-rho_G2)*dAlpha_Lt__dx2;
drho_t__dx3 = (1-Alpha_Lt)*drho_G2__dx3;
drho_t__dx4 = (rho_L-rho_G2)*dAlpha_Lt__dx4 + (1-Alpha_Lt)*drho_G2__dx4;
drho_t__du1 = 0;
drho_t__du2 = (rho_L-rho_G2)*dAlpha_Lt__du2;
drho_t__du3 = (rho_L-rho_G2)*dAlpha_Lt__du3;

%%
Alpha_Lmt = Alpha_Lt*rho_L/(Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2);

dAlpha_Lmt__dx1 = 0;
dAlpha_Lmt__dx2 = (rho_L*dAlpha_Lt__dx2*(Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)...
    - rho_L*Alpha_Lt*dAlpha_Lt__dx2*(rho_L-rho_G2)) ...
    / (Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)^2 ;
dAlpha_Lmt__dx3 = -rho_L*Alpha_Lt*(1-Alpha_Lt)*drho_G2__dx3 /(Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)^2;
dAlpha_Lmt__dx4 = (rho_L*dAlpha_Lt__dx4*(Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)...
    - rho_L*Alpha_Lt*dAlpha_Lt__dx4*(rho_L-rho_G2)) ...
    / (Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)^2 ...
    -rho_L*Alpha_Lt*(1-Alpha_Lt)*drho_G2__dx4 ...
    /(Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)^2;
dAlpha__Lmt__du1 = 0;
dAlpha_Lmt__du2 = (rho_L*dAlpha_Lt__du2*(Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)...
    - rho_L*Alpha_Lt*dAlpha_Lt__du2*(rho_L-rho_G2)) ...
    / (Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)^2 ;
dAlpha_Lmt__du3 = (rho_L*dAlpha_Lt__du3*(Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)...
    - rho_L*Alpha_Lt*dAlpha_Lt__du3*(rho_L-rho_G2)) ...
    / (Alpha_Lt*rho_L + (1-Alpha_Lt)*rho_G2)^2 ;

%%
if P2-P0 <0 
    disp('The flowdirection is wrong. Check parameters 3');
    return
end
w_mix_out = K_pc*u1*sqrt(rho_t*(P2-P0)); 

dw_mix__dx1 = 0;
dw_mix__dx2 = K_pc*u1*(drho_t__dx2*(P2-P0))/(2*sqrt(rho_t*(P2-P0)));
dw_mix__dx3 = K_pc*u1*(drho_t__dx3*(P2-P0)+rho_t*dP2__dx3)/(2*sqrt(rho_t*(P2-P0)));
dw_mix__dx4 = K_pc*u1*(drho_t__dx4*(P2-P0)+rho_t*dP2__dx4)/(2*sqrt(rho_t*(P2-P0)));
dw_mix__du1 = K_pc*sqrt(rho_t*(P2-P0));
dw_mix__du2 = K_pc*u1*(drho_t__du2*(P2-P0))/(2*sqrt(rho_t*(P2-P0)));
dw_mix__du3 = K_pc*u1*(drho_t__du3*(P2-P0))/(2*sqrt(rho_t*(P2-P0)));
%%
% Q_out = K_pc*u1*sqrt((P2-P0)/rho_t); 
% 
% dQ__dx1 = 0;
% dQ__dx2 = K_pc*u1*(-drho_t__dx2*(P2-P0))/(2*rho_t*sqrt(rho_t*(P2-P0)));
% dQ__dx3 = K_pc*u1*(-drho_t__dx3*(P2-P0)+rho_t*dP2__dx3)/(2*rho_t*sqrt(rho_t*(P2-P0)));
% dQ__dx4 = K_pc*u1*(-drho_t__dx4*(P2-P0)+rho_t*dP2__dx4)/(2*rho_t*sqrt(rho_t*(P2-P0)));
% dQ__du1 = K_pc*sqrt((P2-P0)/rho_t);
% dQ__du2 = K_pc*u1*(-drho_t__du2*(P2-P0))/(2*rho_t*sqrt(rho_t*(P2-P0)));
% dQ__du3 = K_pc*u1*(-drho_t__du3*(P2-P0))/(2*rho_t*sqrt(rho_t*(P2-P0)));
%%
w_L_out = Alpha_Lmt*w_mix_out;

dw_L_out__dx1 = 0;
dw_L_out__dx2 = dAlpha_Lmt__dx2*w_mix_out + dw_mix__dx2*Alpha_Lmt;
dw_L_out__dx3 = dAlpha_Lmt__dx3*w_mix_out + dw_mix__dx3*Alpha_Lmt;
dw_L_out__dx4 = dAlpha_Lmt__dx4*w_mix_out + dw_mix__dx4*Alpha_Lmt;
dw_L_out__du1 = dw_mix__du1*Alpha_Lmt;
dw_L_out__du2 = dAlpha_Lmt__du2*w_mix_out + dw_mix__du2*Alpha_Lmt;
dw_L_out__du3 = dAlpha_Lmt__du3*w_mix_out + dw_mix__du3*Alpha_Lmt;
%%
w_G_out = (1 - Alpha_Lmt)*w_mix_out;

dw_G_out__dx1 = 0;
dw_G_out__dx2 = -dAlpha_Lmt__dx2*w_mix_out + dw_mix__dx2*(1-Alpha_Lmt);
dw_G_out__dx3 = -dAlpha_Lmt__dx3*w_mix_out + dw_mix__dx3*(1-Alpha_Lmt);
dw_G_out__dx4 = -dAlpha_Lmt__dx4*w_mix_out + dw_mix__dx4*(1-Alpha_Lmt);
dw_G_out__du1 =  dw_mix__du1*(1-Alpha_Lmt);
dw_G_out__du2 = -dAlpha_Lmt__du2*w_mix_out + dw_mix__du2*(1-Alpha_Lmt);
dw_G_out__du3 = -dAlpha_Lmt__du3*w_mix_out + dw_mix__du3*(1-Alpha_Lmt);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finished calculating the differensials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State equation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
dx1 = u2 - w_G_lp;
dx2 = u3 - w_L_lp;
dx3 = w_G_lp - w_G_out; 
dx4 = w_L_lp -w_L_out;

% jackobian:
%% For A:
dx1__dx1 = -dw_G_lp__dx1;
dx1__dx2 = -dw_G_lp__dx2;
dx1__dx3 = -dw_G_lp__dx3;
dx1__dx4 = -dw_G_lp__dx4;

dx2__dx1 = -dw_L_lp__dx1;
dx2__dx2 = -dw_L_lp__dx2;
dx2__dx3 = -dw_L_lp__dx3;
dx2__dx4 = -dw_L_lp__dx4;

dx3__dx1 = dw_G_lp__dx1 - dw_G_out__dx1;
dx3__dx2 = dw_G_lp__dx2 - dw_G_out__dx2;
dx3__dx3 = dw_G_lp__dx3 - dw_G_out__dx3;
dx3__dx4 = dw_G_lp__dx4 - dw_G_out__dx4;

dx4__dx1 = dw_L_lp__dx1 - dw_L_out__dx1;
dx4__dx2 = dw_L_lp__dx2 - dw_L_out__dx2;
dx4__dx3 = dw_L_lp__dx3 - dw_L_out__dx3;
dx4__dx4 = dw_L_lp__dx4 - dw_L_out__dx4;
%% For B:
dx1__du1 = 0;
dx1__du2 = 1 -dw_G_lp__du2; %1
dx1__du3 = -dw_G_lp__du3;

dx2__du1 = 0;
dx2__du2 = -dw_L_lp__du2;
dx2__du3 = 1 -dw_L_lp__du3; %1

dx3__du1 = dw_G_lp__du1 - dw_G_out__du1;
dx3__du2 = dw_G_lp__du2 - dw_G_out__du2;
dx3__du3 = dw_G_lp__du3 - dw_G_out__du3;

dx4__du1 = dw_L_lp__du1 - dw_L_out__du1;
dx4__du2 = dw_L_lp__du2 - dw_L_out__du2;
dx4__du3 = dw_L_lp__du3 - dw_L_out__du3;
%% For C: 
% Comments:
%y= [ P1 P2 w_mix_out]
%y(1) = P1(x) all derivatives for P1 done before
%y(2) = P2(x) all derivatives for P2 done before
%y(3) = w_mix_out(x,u) all derivatives done before
%% For D
% Comments: 
%only y(3) is a function of u all derivatives done before 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connecting jackobians to form A, B, C and D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [ dx1__dx1 dx1__dx2 dx1__dx3 dx1__dx4 ; 
      dx2__dx1 dx2__dx2 dx2__dx3 dx2__dx4 ;
      dx3__dx1 dx3__dx2 dx3__dx3 dx3__dx4 ; 
      dx4__dx1 dx4__dx2 dx4__dx3 dx4__dx4 ];

B = [ dx1__du1 dx1__du2 dx1__du3 ; 
      dx2__du1 dx2__du2 dx2__du3 ;
      dx3__du1 dx3__du2 dx3__du3 ;
      dx4__du1 dx4__du2 dx4__du3 ];

C = 1e-5.*[ dP1__dx1 dP1__dx2 dP1__dx3 dP1__dx4 ;
            dP2__dx1 dP2__dx2 dP2__dx3 dP2__dx4 ] ; %Pa-> bar
%C(3,:) =  1000*[ dQ__dx1 dQ__dx2 dQ__dx3 dQ__dx4 ]; %[L/s]
C(3,:) =  [ dw_mix__dx1 dw_mix__dx2 dw_mix__dx3 dw_mix__dx4 ]; %[kg/s]

D = zeros(3,3);
D(3,:) = [dw_mix__du1 dw_mix__du2 dw_mix__du3 ];