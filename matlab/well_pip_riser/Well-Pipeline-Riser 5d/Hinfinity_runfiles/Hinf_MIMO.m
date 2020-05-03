function [ GAMMA3 ] = Hinf_MIMO()
%% SISO Topside Valve Z2 %%
cd('..\Hinfinity_MIMO')

%%
run v3_Hinf_MIMO_PbhPin

if isempty(GAM)
    GAMMA3(1) = inf
else
    GAMMA3(1) = GAM
end
save('gamster1.mat','GAMMA3')
%%
run v3_Hinf_MIMO_PbhPrb
load('gamster1.mat')
if isempty(GAM)
    GAMMA3(2) = inf
else
    GAMMA3(2) = GAM
end
save('gamster1.mat','GAMMA3')
%%
run v3_Hinf_MIMO_PbhPt
load('gamster1.mat')
if isempty(GAM)
    GAMMA3(3) = inf
else
    GAMMA3(3) = GAM
end
save('gamster1.mat','GAMMA3')
%%
run v3_Hinf_MIMO_PbhWout
load('gamster1.mat')
if isempty(GAM)
    GAMMA3(4) = inf
else
    GAMMA3(4) = GAM
end
save('gamster1.mat','GAMMA3')
%%
run v3_Hinf_MIMO_PwhPrb
load('gamster1.mat')
if isempty(GAM)
    GAMMA3(5) = inf
else
    GAMMA3(5) = GAM
end
save('gamster1.mat','GAMMA3')

%%
run v3_Hinf_MIMO_PwhPt
load('gamster1.mat')
if isempty(GAM)
    GAMMA3(6) = inf
else
    GAMMA3(6) = GAM
end
save('gamster1.mat','GAMMA3')
%%
run v3_Hinf_MIMO_PwhWout
load('gamster1.mat')
if isempty(GAM)
    GAMMA3(7) = inf
else
    GAMMA3(7) = GAM
end
save('gamster1.mat','GAMMA3')

%%
run v3_Hinf_MIMO_PinPt
load('gamster1.mat')
if isempty(GAM)
    GAMMA3(8) = inf
else
    GAMMA3(8) = GAM
end
save('gamster1.mat','GAMMA3')

end