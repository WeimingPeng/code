%% function:  initialize mode
% university: SWPU 
% time: 2020
% design by WeimingPeng

function [x0, y0] = well_initialize(u0, par)
    %% initialize well&manifold
    ss1 = fsolve(@(x) func(x,par.p_re1,par.glr1,par.rho_l1),[0,0]);
    ss2 = fsolve(@(x) func(x,par.p_re2,par.glr2,par.rho_l2),[0,0]);
    ss3 = fsolve(@(x) func(x,par.p_re3,par.glr3,par.rho_l3),[0,0]);
    ss4 = fsolve(@(x) func(x,par.p_re4,par.glr4,par.rho_l4),[0,0]);
    
    [x1,y1] = initfunc1(ss1,par,u0(1),par.rho_l1,par.glr1,par.p_re1);
    [x2,y2] = initfunc1(ss2,par,u0(2),par.rho_l2,par.glr2,par.p_re2);
    [x3,y3] = initfunc1(ss3,par,u0(3),par.rho_l3,par.glr3,par.p_re3);
    [x4,y4] = initfunc1(ss4,par,u0(4),par.rho_l4,par.glr4,par.p_re4);
    
    y0 = [y1, y2, y3,y4];
    x0 = [x1, x2, x3, x4];
    
    %% function
    function func =func(x, p_re, glr,rho_l)
        rho_m = (x(1)+x(2))/(par.l_w*par.a_w);
        p_wh = par.R*par.t_tw*(x(1)/(par.l_w*par.a_w-x(2)/rho_l))/par.M*1e-5;
        if p_wh>par.p_out
            q_wh = par.c1*u0(1)*sqrt(rho_m*(p_wh-par.p_out));
        else
            q_wh = 0;
        end
        p_bh = p_wh+rho_m*par.l_w*par.g*1e-5;
        q_bhl = par.PI*(p_re-p_bh)*1e5;
        q_bh = q_bhl*(1+glr);
        eq1 = x(1)-glr*x(2);
        eq2 = q_bh-q_wh;

        func = [eq1;eq2];
    end

end
function [x0,y0] = initfunc1(ss,par,u,rho_l,glr,p_re)
    x1_1ss = real(double(ss(1)));
    x1_2ss = real(double(ss(2)));
    rho_mss = (x1_1ss+x1_2ss)/(par.l_w*par.a_w);
    p_whss = par.R*par.t_tw/par.M*x1_1ss/(par.l_w*par.a_w-x1_2ss/rho_l)*1e-5;
    p_bhss = p_whss+rho_mss*par.l_w*par.g*1e-5;
    q_bhss = par.PI*(p_re-p_bhss)*1e5*(1+glr);
    q_whss = par.c1*u*sqrt(rho_mss*(p_whss-par.p_out));
    
    x0(1:2,1) = [x1_1ss; x1_2ss];
    y0(1:4,1) = [p_whss; p_bhss;q_bhss;q_whss];
end

 
