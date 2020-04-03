%%% Save file for Hinf analysis %%%
clc
clear all
here = pwd;
%%
%Set linearization point
z1=0.1; %Subsea valve
z2=0.1; %Topside valve

cd('..\WPR_model')
set_lin_point(z1,z2); % 
%%
%% MIMO Topside Valve %%

%%

cd(here)
Hinf_MIMO_1 = Hinf_MIMO_1()';
%%
xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_MIMO.xls',Hinf_MIMO_1,'Controllability_data_MIMO_1010','I2:I9')


%%
cd(here)
Hinf_MIMO_2 = Hinf_MIMO_2()';
%%
xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_MIMO.xls',Hinf_MIMO_2,'Controllability_data_MIMO_1010','J2:J9')
%%

cd(here)
Hinf_MIMO = Hinf_MIMO()';


%%
xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_MIMO.xls',Hinf_MIMO,'Controllability_data_MIMO_1010','K2:K9')