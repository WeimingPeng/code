%% function:  initialize sulg
% university: SWPU
% time: 2020
% design by WeimingPeng

%{ 
use nominal mass to initialize pip_riser_top
%}

%% main code
function [x1, y1] = a_initialize_slug(u, par)
    z=u(1);
    wG_in=u(2);
    wL_in=u(3);
    
    syms P1 x3 x4
    h1 = 0.07;
    P2_t = x3*par.R*par.T2/(par.M_G*(par.V2 - x4/par.rho_L));
    rho_G2 = P2_t*par.M_G/(par.R*par.T2);
    Fric_pipe = 1.0164e+005;
    Fric_riser = 3.6497e+004;
    Alpha_L2_av = x4/(par.rho_L*par.V2);

    rho_mix_av = (x3+x4)/par.V2;
    P2_b = P2_t + rho_mix_av*par.g*par.L2 + Fric_riser;
    A_g = par.A1 + par.A1*h1^2/par.hc^2 - 2*h1*par.A1/par.hc;
    A_l = par.A1-A_g ;

    eq1 = wL_in^2 - par.rho_L*(P1-Fric_pipe+par.rho_L*par.g*h1-P2_b)*(A_l*par.K_o)^2;

    Alpha_Lt =  2*Alpha_L2_av - h1/par.hc;
    rho_t = Alpha_Lt*par.rho_L + (1-Alpha_Lt)*rho_G2;
    eq2 = (wL_in+wG_in)^2 - rho_t*(P2_t-par.P0)*(par.K_pc*z)^2;
    eq3 = wL_in/(wG_in + wL_in) - x4/(x4+x3);
    %%
    ss = solve(eq1,eq2,eq3);
    %%
    if(strcmp('7.4.0.287 (R2007a)',version)) %tf = strcmp(s1,s2) 比较 s1 和 s2，如果二者相同，则返回 1 (true)，否则返回 0 (false)。
        ind = 3;
    else
        ind = 2;
    end
    P1_ss = real(double(ss.P1(ind))); % real返回数组中的实部
    x3_ss = real(double(ss.x3(ind)));
    x4_ss = real(double(ss.x4(ind)));
    
    rho_G1_n = P1_ss*par.M_G/(par.R*par.T1);
    a_new = wL_in*rho_G1_n/(wL_in*rho_G1_n + wG_in*par.rho_L);% alpha_lm即，体积分数
    a_old = a_new - 0.01; % to make sure loop runs at least once

    h1_old = h1;
    h1_new = par.k_h*a_new*par.hc; %mean_h1

    opt = optimset('Display','off','TolFun',1e-20,'TolX',1e-20,'MaxIter',1e4,'MaxFunEvals',1e10);
    y0(1) = P1_ss/1e5; % 单位转换
    while(abs(a_new-a_old)>1e-13)
        par.P1_norm = 1e5*y0(1);
        rho_G1 = 1e5*y0(1)*par.M_G/(par.R*par.T1);
        par.Alpha_L1_av = wL_in*rho_G1/(wL_in*rho_G1 + wG_in*par.rho_L);
        x2_ss = par.rho_L*par.V1*a_new + (h1_old - h1_new)*par.A1*par.rho_L*(1-a_new)/sin(par.theta);
        x1_ss = (par.V1 - x2_ss/par.rho_L)*P1_ss*par.M_G/(par.R*par.T1);
        x0 = [x1_ss;x2_ss;x3_ss;x4_ss];
        x0 = fsolve(@a_model0_slug,x0,opt,u,par);
        y0 = a_model_slug(0,x0,u,'measurements',par);
        rho_G1 = 1e5*y0(1)*par.M_G/(par.R*par.T1);
        a_old = a_new;
        a_new = wL_in*rho_G1/(wL_in*rho_G1 + wG_in*par.rho_L);
        h1_old = h1_new;
        h1_new = par.k_h*a_new*par.hc;
    end
 
    x1 = x0;
    y1 = y0;

    str1 = ['Pressure in pipeline: ' num2str(y0(1)) ' bar' ];
    disp(str1)
    str2 = ['Pressure at top of riser: ' num2str(y0(2)) ' bar' ];
    disp(str2)
    str2 = ['Flow rate of choke: ' num2str(y0(3)) ' kg/s' ];
    disp(str2)
end