function [ GAMMA1 ] = Hinf_SISO_Z1()
%% SISO Subsea Valve Z1 %%
cd('..\Hinfinity_Z1')

%%
run v3_Hinf_P_bh_Z1

if isempty(GAM)
    GAMMA1(1) = inf
else
    GAMMA1(1) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_P_wh_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(2) = inf
else
    GAMMA1(2) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_W_in_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(3) = inf
else
    GAMMA1(3) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_P_in_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(4) = inf
else
    GAMMA1(4) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_P_rb_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(5) = inf
else
    GAMMA1(5) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_DPr_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(6) = inf
else
    GAMMA1(6) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_P_t_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(7) = inf
else
    GAMMA1(7) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_Q_out_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(8) = inf
else
    GAMMA1(8) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_W_out_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(9) = inf
else
    GAMMA1(9) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_Rho_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(10) = inf
else
    GAMMA1(10) = GAM
end
save('gamster1.mat','GAMMA1')
%%
run v3_Hinf_Alpha_Z1
load('gamster1.mat')
if isempty(GAM)
    GAMMA1(11) = inf
else
    GAMMA1(11) = GAM
end
save('gamster1.mat','GAMMA1')


end

