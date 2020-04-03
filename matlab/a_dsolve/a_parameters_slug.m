%% function:  parameter slug
% university: SWPU
% time: 2020
% design by WeimingPeng


%% main code
function par = a_parameters_slug()
    
    par.g = 9.81;                       % Gravity (m/s2)
    par.R = 8314;                       % Gas constant (J/(K.Kmol))
    par.rho_L = 805;                  % oil density (Kg/m3)

    par.theta= pi/180;                  % Feed pipe inclination (Rad) = 1 Deg
    par.r1 = 0.06;                      % Raduis of pipe (m)
    par.r2 = 0.05;                      % Raduis of riser (m)
    par.hc = 2*par.r1/cos(par.theta);   % Critical liquid level (m)
    par.A1 = pi*par.r1^2;               % Cross section area of pipeline (m2)
    par.A2 = pi*par.r2^2;               % Cross section area of riser (m2)
    par.L1 = 4300;                      % Length of upstream pipe (m)
    par.L2 = 300;                       % Height of riser
    par.L3 = 100;                       % length of horozontal top section (m)
    par.T1 = 335;                       % Pipeline temprature (K)
    par.T2 = 298.3;                     % Riser temprature (K)
    par.M_G = 23;                    % Molecular weight of Gas (kg/kmol)
    par.P0 = 50e5;                    % Pressure after choke valve (Pa)
    par.V2 = par.A2*(par.L2+par.L3);    % Total volume in riser (m3)
    par.V1 = par.A1*par.L1;             % Total volume in upstream (m3)
    %% *************************************************************************
    par.Pnormal = 101325;              % atmospheric conditions for Usg
    Tnormal = 293.15;                  % Temprature in Normal conditions
    RTnormal =par.R*Tnormal/par.M_G;   % compressibility in Normal Conditions
    par.Rog_n = par.Pnormal/RTnormal;  % Gas density in normal conditions
    par.visl= 1.4260e-4;
    par.eps = 2.8e-5;                  % Roughness of pipe

    %% ********************************************************
    %% constant // 使用OLGA计算结果计算几个常数得值
    z0= 0.04;
    P1 = 76.59951e5;
    P2_t = 57.81898e5;
    wG_in = 0.36;
    wL_in = 8.64;
    
    w_G1 = wG_in;
    w_L1 = wL_in;
    wG_out = wG_in;
    wL_out = wL_in;
    w_mix = wG_out + wL_out;

    rho_G1 = P1*par.M_G/(par.R*par.T1);
    Alpha_L1_av = wL_in*rho_G1/(wL_in*rho_G1 + wG_in*par.rho_L);

    par.k_h = 0.74;
    h1ss = par.k_h*Alpha_L1_av*par.hc;
    h1 = h1ss;
    Uslin = wL_in/(par.A1*par.rho_L);

    Re1=par.rho_L*Uslin*(2*par.r1)/par.visl;
    Lambda1=0.0056 + 0.5*Re1^(-0.32);

    Fric_pipe=0.5*Alpha_L1_av*Lambda1*par.rho_L*Uslin^2*par.L1/(2*par.r1);

    P2_av = (P1+par.rho_L*par.g*h1-Fric_pipe+P2_t)/2;
    rho_G2_av = (P2_av)*par.M_G/(par.R*par.T2);
    Alpha_L2_av = rho_G2_av*wL_in/(rho_G2_av*wL_in+par.rho_L*wG_in);
    rho_mix_av = rho_G2_av*(1-Alpha_L2_av)+par.rho_L*Alpha_L2_av;

    Usl2 = wL_in/(par.rho_L*par.A2);
    Usg2 = wG_in/(rho_G2_av*par.A2);
    Um = Usl2 + Usg2;
    Re2=rho_mix_av*Um*2*par.r2/par.visl;
    Lambda2=0.0056 + 0.5*Re2^(-0.32);
    Fric_riser=0.5*Alpha_L2_av*Lambda2*rho_mix_av*Um^2*(par.L2+par.L3)/(2*par.r2);

    P2_b = P2_t + rho_mix_av*par.g*par.L2 + Fric_riser;

    rho_G2 = P2_t*par.M_G/(par.R*par.T2);
    Alpha_Lt = rho_G2*w_L1/(rho_G2*w_L1+par.rho_L*w_G1);

    rho_t = Alpha_Lt*par.rho_L + (1-Alpha_Lt)*rho_G2;

    A_g = (h1<par.hc)*(par.A1/par.hc^2)*(par.hc - h1)^2;
    A_l = par.A1 - A_g;

    v_G1 = w_G1/(rho_G1*A_g);
    par.K_g = 1.207*v_G1/sqrt((P1-Fric_pipe-P2_b)/rho_G1);
    par.K_o = 1.1*w_L1/(A_l*sqrt(par.rho_L*(P1-Fric_pipe+par.rho_L*par.g*h1-P2_b)));
    par.K_pc = 1*w_mix/(z0*sqrt(rho_t*(P2_t-par.P0)));
end