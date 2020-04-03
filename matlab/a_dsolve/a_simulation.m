%% function:  simulation
% university: SWPU
% time: 2020
% design by WeimingPeng

%{ 
1.you can  chose the nominal pressure in separater or nominal mass flow from
reservoir to initianal the mode!
%}

clc
clear
%% main code
tic
disp("initialize well&manifold&separator")
u0 = [0.06, 0.06, 0.06, 0.06, 0.05]; % well&top valve opening
par = a_parameter();
disp("###########")
[x0, y0] = well_initialize(u0, par);
w_gas = (y0(4,1)*par.x1(1)+y0(4,2)*par.x2(1)+y0(4,3)*par.x3(1)+y0(4,4)*par.x4(1));
w_water = (y0(4,1)*par.x1(2)+y0(4,2)*par.x2(2)+y0(4,3)*par.x3(2)+y0(4,4)*par.x4(2));
w_oil = (y0(4,1)*par.x1(3)+y0(4,2)*par.x2(3)+y0(4,3)*par.x3(3)+y0(4,4)*par.x4(3));
w_sep = y0(4,1)+y0(4,2)+y0(4,3)+y0(4,4);
gas_ratio = w_gas/w_sep;
water_ratio = w_water/w_sep;
oil_ratio = w_oil/w_sep;
glr = gas_ratio/(water_ratio+oil_ratio);
wc = water_ratio/(water_ratio+oil_ratio);
rho_l = par.rho_oil*par.rho_water*(oil_ratio+water_ratio)/(par.rho_oil*water_ratio+par.rho_water*oil_ratio);

z0 = separator_initialize(par,glr,wc,rho_l);
w_nom(1) = w_sep*gas_ratio;
w_nom(2) = w_sep*oil_ratio;

u_slug = [u0(5),w_nom(1),w_nom(2)];
par1 = a_parameters_slug();
[x1, y1] = a_initialize_slug(u_slug, par1);
% 
% % initialize pump station calculate \Delta P
delta_p = 30; % (bar)
p_pumpb = y1(1)-delta_p;
dpg = (par.p_out-p_pumpb)*1e5;
rho_sepg = par.p_out*par.M*1e5/par.R/par.t_tw;
dpl = par.p_out*1e5+par.rho_oil*par.g*(z0(4)+z0(5))-p_pumpb*1e5;
z1_g = w_nom(1)/sqrt(dpg*rho_sepg)/par.cv_sepg;
z2_l = w_nom(2)/sqrt(dpl*par.rho_oil)/par.cv_sepl; % 阀门公式没有用分离器的公式，貌似哪里有问题
% 
% % input
% u = {u0(1:4);[z1_g;z2_l;delta_p];u0(5)};
% % state
% % x = 
% % output initial
% % y = 
% 
disp("initialize complete")
toc

