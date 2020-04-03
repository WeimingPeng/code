function [GAMMA_SIMO1] = Hinf_SIMO_Z1()
%% SIMO Subsea Valve Z1 %%
cd('..\Hinfinity_Z1')

%%
run v3_Hinf_SIMO_Z1_PbhPwh

if isempty(GAM)
    GAMMA_SIMO1(1) = inf
else
    GAMMA_SIMO1(1) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO1')
%%
run v3_Hinf_SIMO_Z1_PbhWin
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO1(2) = inf
else
    GAMMA_SIMO1(2) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO1')
%%
run v3_Hinf_SIMO_Z1_PbhPin
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO1(3) = inf
else
    GAMMA_SIMO1(3) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO1')
%%
run v3_Hinf_SIMO_Z1_PwhWin
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO1(4) = inf
else
    GAMMA_SIMO1(4) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO1')
%%
run v3_Hinf_SIMO_Z1_PwhPin
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO1(5) = inf
else
    GAMMA_SIMO1(5) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO1')
%%
run v3_Hinf_SIMO_Z1_WinPin
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO1(6) = inf
else
    GAMMA_SIMO1(6) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO1')

end