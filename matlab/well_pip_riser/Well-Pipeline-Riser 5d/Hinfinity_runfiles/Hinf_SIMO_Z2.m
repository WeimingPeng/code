function [ GAMMA_SIMO2 ] = Hinf_SIMO_Z2(  )
%% SIMO Topside Valve Z2 %%
cd('..\Hinfinity_Z2')

%%
run v3_Hinf_SIMO_Z2_PinPt

if isempty(GAM)
    GAMMA_SIMO2(1) = inf
else
    GAMMA_SIMO2(1) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO2')
%%
run v3_Hinf_SIMO_Z2_PinWout
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO2(2) = inf
else
    GAMMA_SIMO2(2) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO2')
%%
run v3_Hinf_SIMO_Z2_PtWout
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO2(3) = inf
else
    GAMMA_SIMO2(3) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO2')
%%
run v3_Hinf_SIMO_Z2_PtRho
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO2(4) = inf
else
    GAMMA_SIMO2(4) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO2')
%%
run v3_Hinf_SIMO_Z2_PbhWout
load('gamsterSIMO.mat')
if isempty(GAM)
    GAMMA_SIMO2(5) = inf
else
    GAMMA_SIMO2(5) = GAM
end
save('gamsterSIMO.mat','GAMMA_SIMO2')


end

