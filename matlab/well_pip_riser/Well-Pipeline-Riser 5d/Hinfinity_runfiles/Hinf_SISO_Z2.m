function [ GAMMA2 ] = Hinf_SISO_Z2()
%% SISO Topside Valve Z2 %%
cd('..\Hinfinity_Z2')

%%
run v3_Hinf_P_bh_Z2

if isempty(GAM)
    GAMMA2(1) = inf
else
    GAMMA2(1) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_P_wh_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(2) = inf
else
    GAMMA2(2) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_W_in_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(3) = inf
else
    GAMMA2(3) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_P_in_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(4) = inf
else
    GAMMA2(4) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_P_rb_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(5) = inf
else
    GAMMA2(5) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_DPr_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(6) = inf
else
    GAMMA2(6) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_P_t_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(7) = inf
else
    GAMMA2(7) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_Q_out_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(8) = inf
else
    GAMMA2(8) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_W_out_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(9) = inf
else
    GAMMA2(9) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_Rho_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(10) = inf
else
    GAMMA2(10) = GAM
end
save('gamster1.mat','GAMMA2')
%%
run v3_Hinf_Alpha_Z2
load('gamster1.mat')
if isempty(GAM)
    GAMMA2(11) = inf
else
    GAMMA2(11) = GAM
end
save('gamster1.mat','GAMMA2')


end

