

clc
clear
par = a_parameter();
rho_gsep = par.p_out*par.M*1e5/par.R/par.t_tw;
A = [0.036*par.rho_l,-rho_gsep;1,1];
B = [0;pi*par.r_sep^2*par.l_sep];
x = linsolve(A,B);
x
xxx=x(2)*rho_gsep/(par.rho_l*x(1))

    function func1 = func1(xse)
        rho_gsep = par.p_out*par.M*1e5/par.R/par.t_tw;
        eq11 = rho_gsep*xse(1)/par.rho_l-par.glr1;
        eq12 = xse(1)+xse(2)-pi*par.r_sep^2*par.l_sep;
        func1 = [eq11;eq12];
    end


   