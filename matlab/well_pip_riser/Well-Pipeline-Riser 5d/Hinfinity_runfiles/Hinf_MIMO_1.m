function [ GAMMA4 ] = Hinf_MIMO_1()
%% SISO Topside Valve Z2 %%
cd('..\Hinfinity_MIMO')

%%
run v3_Hinf_MIMO_PbhPin_1

if isempty(GAM)
    GAMMA4(1) = inf
else
    GAMMA4(1) = GAM
end
save('gamster2.mat','GAMMA4')
%%
run v3_Hinf_MIMO_PbhPrb_1
load('gamster2.mat')
if isempty(GAM)
    GAMMA4(2) = inf
else
    GAMMA4(2) = GAM
end
save('gamster2.mat','GAMMA4')
%%
run v3_Hinf_MIMO_PbhPt_1
load('gamster2.mat')
if isempty(GAM)
    GAMMA4(3) = inf
else
    GAMMA4(3) = GAM
end
save('gamster2.mat','GAMMA4')
%%
run v3_Hinf_MIMO_PbhWout_1
load('gamster2.mat')
if isempty(GAM)
    GAMMA4(4) = inf
else
    GAMMA4(4) = GAM
end
save('gamster2.mat','GAMMA4')
%%
run v3_Hinf_MIMO_PwhPrb_1
load('gamster2.mat')
if isempty(GAM)
    GAMMA4(5) = inf
else
    GAMMA4(5) = GAM
end
save('gamster2.mat','GAMMA4')

%%
run v3_Hinf_MIMO_PwhPt_1
load('gamster2.mat')
if isempty(GAM)
    GAMMA4(6) = inf
else
    GAMMA4(6) = GAM
end
save('gamster2.mat','GAMMA4')
%%
run v3_Hinf_MIMO_PwhWout_1
load('gamster2.mat')
if isempty(GAM)
    GAMMA4(7) = inf
else
    GAMMA4(7) = GAM
end
save('gamster2.mat','GAMMA4')

%%
run v3_Hinf_MIMO_PinPt_1
load('gamster2.mat')
if isempty(GAM)
    GAMMA4(8) = inf
else
    GAMMA4(8) = GAM
end
save('gamster2.mat','GAMMA4')

end