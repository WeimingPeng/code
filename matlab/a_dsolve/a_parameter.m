%% function:  initial parameter
% university: SWPU
% time: 2020
% design by WeimingPeng
% 

function par = a_parameter( )
    %% fundamental parameter
    par.p_out = 60; % bar p_se
    par.g = 9.81;
    
    %% reservior and wells parameter
    par.p_re1 = 282; % well1 pre(bar)
    par.p_re2 = 260;
    par.p_re3 = 270;
    par.p_re4 = 265;
    par.wc1= 0.7; % well1 water cut
    par.wc2= 0.6; 
    par.wc3= 0.6; 
    par.wc4= 0.4; 
    
    par.gor1 = 0.09;
    par.gor2 = 0.08;
    par.gor3 = 0.09;
    par.gor4 = 0.07;

   par.d_w = 0.12; % well diameter(m)
   par.a_w = pi*par.d_w^2/4;
   par.l_w = 2000; % Height of oil well (m)
   par.t_tw= 350; % Temprature well head(top well,K)
   par.PI = 1.34e-6; % Productivity Index of the well
   par.M = 0.02006; %(kg/mol)
   par.R = 8.3145;
   par.c1 = 0.506; % k1 valve coefficient
   par.rho_oil = 900; % oil density (Kg/m3)
   par.rho_water = 1000; % water desity(Kg/m3)
   [par.glr1,par.rho_l1,x1] = parfunc1(par.gor1, par.wc1, par.rho_oil, par.rho_water);
   [par.glr2,par.rho_l2,x2] = parfunc1(par.gor2, par.wc2, par.rho_oil, par.rho_water);
   [par.glr3,par.rho_l3,x3] = parfunc1(par.gor3, par.wc3, par.rho_oil, par.rho_water);
   [par.glr4,par.rho_l4,x4] = parfunc1(par.gor4, par.wc4, par.rho_oil, par.rho_water);

   
   %% separator parameter
   par.x1 = x1;
   par.x2 = x2;
   par.x3 = x3;
   par.x4 = x4;  
%    par.beta_gr = mean([x1(1),x2(1),x3(1),x4(1)]);
%    par.beta_wr = mean([x1(2),x2(2),x3(2),x4(2)]);
%    par.beta_or = mean([x1(3),x2(3),x3(3),x4(3)]); % gas,oil,water ratio
   par.r_sep = 2; %(m)
   par.l_sep = 9; %(m)
   par.h_weir = 2.2; %(m)
   par.l_weir = 2; %(m)
   par.l_water = par.l_sep-par.l_weir;
   par.cv_sepg = 0.0014;  % gas valve constant of separator 
   par.cv_sepl = 0.0014; %% oil valve constant of separator 

end

%% subfunction
    function [glr,rho_l,x] =parfunc1(gor, wc, rho_oil, rho_water)
       B1 = [0;0;1];
       A1 = [1,0,-gor;0,1-wc,-wc;1,1,1];
       x = linsolve(A1,B1);
       glr = x(1)/(x(2)+x(3)); % mass ratio
       rho_l = rho_oil*rho_water*(x(2)+x(3))/(rho_oil*x(2)+rho_water*x(3));  %liquid density(kg/m3)
    end
