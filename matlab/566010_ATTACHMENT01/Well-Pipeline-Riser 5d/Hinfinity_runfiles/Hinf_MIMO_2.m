function [ GAMMA5 ] = Hinf_MIMO_2()
%% SISO Topside Valve Z2 %%
cd('..\Hinfinity_MIMO')

%%
run v3_Hinf_MIMO_PbhPin_2

if isempty(GAM)
    GAMMA5(1) = inf
else
    GAMMA5(1) = GAM
end
save('gamster3.mat','GAMMA5')
%%
run v3_Hinf_MIMO_PbhPrb_2
load('gamster3.mat')
if isempty(GAM)
    GAMMA5(2) = inf
else
    GAMMA5(2) = GAM
end
save('gamster3.mat','GAMMA5')
%%
run v3_Hinf_MIMO_PbhPt_2
load('gamster3.mat')
if isempty(GAM)
    GAMMA5(3) = inf
else
    GAMMA5(3) = GAM
end
save('gamster3.mat','GAMMA5')
%%
run v3_Hinf_MIMO_PbhWout_2
load('gamster3.mat')
if isempty(GAM)
    GAMMA5(4) = inf
else
    GAMMA5(4) = GAM
end
save('gamster3.mat','GAMMA5')
%%
run v3_Hinf_MIMO_PwhPrb_2
load('gamster3.mat')
if isempty(GAM)
    GAMMA5(5) = inf
else
    GAMMA5(5) = GAM
end
save('gamster3.mat','GAMMA5')

%%
run v3_Hinf_MIMO_PwhPt_2
load('gamster3.mat')
if isempty(GAM)
    GAMMA5(6) = inf
else
    GAMMA5(6) = GAM
end
save('gamster3.mat','GAMMA5')
%%
run v3_Hinf_MIMO_PwhWout_2
load('gamster3.mat')
if isempty(GAM)
    GAMMA5(7) = inf
else
    GAMMA5(7) = GAM
end
save('gamster3.mat','GAMMA5')

%%
run v3_Hinf_MIMO_PinPt_2
load('gamster3.mat')
if isempty(GAM)
    GAMMA5(8) = inf
else
    GAMMA5(8) = GAM
end
save('gamster3.mat','GAMMA5')

end