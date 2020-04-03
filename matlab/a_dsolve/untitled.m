

ss = fsolve(@func_slug);
x1 = ss;
y1 = ss;
function func_slug= func_slug(x)
    par = a_parameters_slug();
    z=0.04;
    w_nom = [0.3,5.84];
    wG_in=w_nom(1);
    wL_in=w_nom(2);

    h1 = 0.07;
    P2_t = x(1)*par.R*par.T2/(par.M_G*(par.V2 - x(2)/par.rho_L));
    rho_G2 = P2_t*par.M_G/(par.R*par.T2);
    Fric_pipe = 1.0164e+005;
    Fric_riser = 3.6497e+004;
    Alpha_L2_av = x(2)/(par.rho_L*par.V2);

    rho_mix_av = (x(1)+x(2))/par.V2;
    P2_b = P2_t + rho_mix_av*par.g*par.L2 + Fric_riser;
    A_g = par.A1 + par.A1*h1^2/par.hc^2 - 2*h1*par.A1/par.hc;
    A_l = par.A1-A_g ;

    eq1 = wL_in^2 - par.rho_L*(x(3)-Fric_pipe+par.rho_L*par.g*h1-P2_b)*(A_l*par.K_o)^2;

    Alpha_Lt =  2*Alpha_L2_av - h1/par.hc;
    rho_t = Alpha_Lt*par.rho_L + (1-Alpha_Lt)*rho_G2;
    eq2 = (wL_in+wG_in)^2 - rho_t*(P2_t-par.P0)*(par.K_pc*z)^2;
    eq3 = wL_in/(wG_in + wL_in) - x(2)/(x(2)+x(1));
    
    func_slug = [eq1;eq2;eq3];
end