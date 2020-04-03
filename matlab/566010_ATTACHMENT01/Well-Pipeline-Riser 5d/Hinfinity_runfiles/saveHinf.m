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

%% SISO Subsea Valve Z1 %%
cd(here)
Hinf_SISO_Z1 = Hinf_SISO_Z1()'



% SIMO Subsea Valve Z1 %%
cd(here)
Hinf_SIMO_Z1 = Hinf_SIMO_Z1()'
%

% SISO Topside Valve Z2 %%
cd(here)
Hinf_SISO_Z2 = Hinf_SISO_Z2()'
%

% SIMO Topside Valve Z2 %%
cd(here)
Hinf_SIMO_Z2 = Hinf_SIMO_Z2()'
%%
cd(here)
xlswrite('Controllability tables\Controllability_data_Z1_1.xls',Hinf_SISO_Z1,'Controllability_data_Z1_1010','M2:M12')
xlswrite('Controllability tables\Controllability_data_Z1_1.xls',Hinf_SIMO_Z1,'Controllability_data_Z1_1010','M13:M18')
xlswrite('Controllability tables\Controllability_data_Z2_1.xls',Hinf_SISO_Z2,'Controllability_data_Z2_1010','M2:M12')
xlswrite('Controllability tables\Controllability_data_Z2_1.xls',Hinf_SIMO_Z2,'Controllability_data_Z2_1010','M13:M17')
%% MIMO Topside Valve %%
%not done


%xlswrite('Controllability_data_Z1_1.xls',contr_var,'Controllability_data_Z1_3010','XXX') % Save SISO Hinf
% 
% save('Hinf_data.mat','Hinf_SISO_Z1');
% 
% load('Hinf_data.mat')
% save('Hinf_data.mat','Hinf_SIMO_Z1', '-append');

% C:\Users\matslie\Dropbox\Vår 2012\Matlab\Well-Pipeline-Riser_mats\Controllability_Z2