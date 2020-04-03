function z0 = separator_initialize(par,glr,wc,rho_l)
  %% initialize separator
   
    rho_gsep = par.p_out*par.M*1e5/par.R/par.t_tw;
    A_sep = [glr*rho_l,-rho_gsep;1,1];
    B_sep = [0;pi*par.r_sep^2*par.l_sep];
    ss1 = linsolve(A_sep,B_sep);
    
    v_g_sep = ss1(2); % gas volume in separator
    v_w_sep = ss1(1)*rho_l*wc/par.rho_water;
    v_o_sep = ss1(1)-v_w_sep;
    rho_gsep1 = par.p_out*par.M*1e5/par.R/par.t_tw;
    x2_1 = rho_gsep1*v_g_sep;
    x2_2 = par.rho_water*v_w_sep;
    x2_3 = par.rho_oil*v_o_sep;
    z0(1:3,1) = [x2_1;x2_2;x2_3];
    
    % calculate h_water
    h_water = func2(v_w_sep,par.l_water);
    z0(4,1) = h_water;
    
      % calculate h_water
    v_l_sep = v_w_sep+v_o_sep;
    h_liquid = func2(v_l_sep,par.l_sep);
    z0(5,1) = h_liquid -h_water; %% h_oil
end
function func2 = func2(v_w_sep,l)
        syms h
        a = v_w_sep/l;
        y = sqrt(2*2*h-h^2);
        c = 2*int(y,0,h);
        eq_sep1 = c-a;
        yy = solve(eq_sep1,h);
        h11 = real(yy);
        func2 = h11;
    end