%%% Save file for controllabillity analysis %%%

z1 = 0.3;
z2 = 0.1;
cd('WPR_model')
set_lin_point(z1,z2); % 

%% SISO Subsea Valve Z1 %%
D_y = spdiags(De)'; % Scaling values

contr_var = [y0(1:11)'; % Values
            D_y; % Scaling
            freqresp(G1,0),freqresp(G2,0),freqresp(G3,0),freqresp(G4,0),freqresp(G5,0),freqresp(G6,0),freqresp(G7,0),freqresp(G8,0),freqresp(G9,0),freqresp(G10,0),freqresp(G11,0); % Steady-state gains
            min(abs(Yp1)), min(abs(Yp2)), min(abs(Yp3)), min(abs(Yp4)), min(abs(Yp5)), min(abs(Yp6)), min(abs(Yp7)), min(abs(Yp8)), min(abs(Yp9)), min(abs(Yp10)), min(abs(Yp11));  % Pole vectors
            Msmin1, Msmin2, Msmin3, Msmin4, Msmin5, Msmin6, Msmin7, Msmin8, Msmin9, Msmin10, Msmin11; % Bounds on S and T
            KSmin1, KSmin2, KSmin3, KSmin4, KSmin5, KSmin6, KSmin7, KSmin8, KSmin9, KSmin10, KSmin11; % Bounds on KS
            SGmin1, SGmin2, SGmin3, SGmin4, SGmin5, SGmin6, SGmin7, SGmin8, SGmin9, SGmin10, SGmin11; % Bounds on SG
            KSGd11, KSGd21, KSGd31, KSGd41, KSGd51, KSGd61, KSGd71, KSGd81, KSGd91, KSGd101, KSGd111; % Bounds on KSGd1
            KSGd12, KSGd22, KSGd32, KSGd42, KSGd52, KSGd62, KSGd72, KSGd82, KSGd92, KSGd102, KSGd112; % Bounds on KSGd2
            SGdmin11, SGdmin21, SGdmin31, SGdmin41, SGdmin51, SGdmin61, SGdmin71, SGdmin81, SGdmin91, SGdmin101, SGdmin111; % Bounds on SGd1
            SGdmin12, SGdmin22, SGdmin32, SGdmin42, SGdmin52, SGdmin62, SGdmin72, SGdmin82, SGdmin92, SGdmin102, SGdmin112 % Bounds on SGd2
]';



xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_Z1_1.xls',contr_var,'Controllability_data_Z1_1010','B2') % Save SISO controllability analysis

%% SIMO Subsea Valve Z1 %%
contr_var_SIMO_Z1 = [Msmin12, Msmin13, Msmin14, Msmin23, Msmin24, Msmin34;                                               % bounds on S and T
                    KSmin12, KSmin13, KSmin14, KSmin23, KSmin24, KSmin34;                                                % Bounds on KS
                    GammaSG12_min, GammaSG13_min, GammaSG14_min, GammaSG23_min, GammaSG24_min, GammaSG34_min;            % Bounds on SG
                    KSGdmin121, KSGdmin131, KSGdmin141, KSGdmin231, KSGdmin241, KSGdmin341;                              % Bounds on KSGd1
                    KSGdmin122, KSGdmin132, KSGdmin142, KSGdmin232, KSGdmin242, KSGdmin342;                              % Bounds on KSGd2 (tight)
                    GammaSGd121_min, GammaSGd131_min, GammaSGd141_min, GammaSGd231_min, GammaSGd241_min, GammaSGd341_min;% Bounds on SGd1
                    GammaSGd122_min, GammaSGd132_min, GammaSGd142_min, GammaSGd232_min, GammaSGd242_min, GammaSGd342_min % Bounds on SGd2
    ]';
    


xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_Z1_1.xls',contr_var_SIMO_Z1,'Controllability_data_Z1_1010','F13') % Save SISO controllability analysis


%% SISO Topside Valve Z2 %%
D_y = spdiags(De)'; % Scaling values

contr_var = [y0(1:11)'; % Values
            D_y; % Scaling
            freqresp(G1,0),freqresp(G2,0),freqresp(G3,0),freqresp(G4,0),freqresp(G5,0),freqresp(G6,0),freqresp(G7,0),freqresp(G8,0),freqresp(G9,0),freqresp(G10,0),freqresp(G11,0); % Steady-state gains
            min(abs(Yp1)), min(abs(Yp2)), min(abs(Yp3)), min(abs(Yp4)), min(abs(Yp5)), min(abs(Yp6)), min(abs(Yp7)), min(abs(Yp8)), min(abs(Yp9)), min(abs(Yp10)), min(abs(Yp11));  % Pole vectors
            Msmin1, Msmin2, Msmin3, Msmin4, Msmin5, Msmin6, Msmin7, Msmin8, Msmin9, Msmin10, Msmin11; % Bounds on S and T
            KSmin1, KSmin2, KSmin3, KSmin4, KSmin5, KSmin6, KSmin7, KSmin8, KSmin9, KSmin10, KSmin11; % Bounds on KS
            SGmin1, SGmin2, SGmin3, SGmin4, SGmin5, SGmin6, SGmin7, SGmin8, SGmin9, SGmin10, SGmin11; % Bounds on SG
            KSGd11, KSGd21, KSGd31, KSGd41, KSGd51, KSGd61, KSGd71, KSGd81, KSGd91, KSGd101, KSGd111; % Bounds on KSGd1
            KSGd12, KSGd22, KSGd32, KSGd42, KSGd52, KSGd62, KSGd72, KSGd82, KSGd92, KSGd102, KSGd112; % Bounds on KSGd2
            SGdmin11, SGdmin21, SGdmin31, SGdmin41, SGdmin51, SGdmin61, SGdmin71, SGdmin81, SGdmin91, SGdmin101, SGdmin111; % Bounds on SGd1
            SGdmin12, SGdmin22, SGdmin32, SGdmin42, SGdmin52, SGdmin62, SGdmin72, SGdmin82, SGdmin92, SGdmin102, SGdmin112 % Bounds on SGd2
]';


xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_Z2_1.xls',contr_var,'Controllability_data_Z2_1010','B2') % Save SISO controllability analysis

%% SIMO Topside Valve Z2 %%
contr_var_SIMO_Z2 = [Msmin47, Msmin19, Msmin49, Msmin79, Msmin710;                                       % bounds on S and T
                    KSmin47, KSmin19, KSmin49, KSmin79, KSmin710;                                        % Bounds on KS
                    GammaSG47_min, GammaSG19_min, GammaSG49_min, GammaSG79_min, GammaSG710_min;          % Bounds on SG
                    KSGdmin471, KSGdmin191, KSGdmin491, KSGdmin791, KSGdmin7101;                         % Bounds on KSGd1
                    KSGdmin472, KSGdmin192, KSGdmin492, KSGdmin792, KSGdmin7102;                         % Bounds on KSGd2 (tight)
                    GammaSGd471_min, GammaSGd191_min, GammaSGd491_min, GammaSGd791_min, GammaSGd7101_min;% Bounds on SGd1
                    GammaSGd472_min, GammaSGd192_min, GammaSGd492_min, GammaSGd792_min, GammaSGd7102_min % Bounds on SGd2
    ]';

xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_Z2_1.xls',contr_var_SIMO_Z2,'Controllability_data_Z2_1010','F13') % Save SISO controllability analysis


%% MIMO Topside Valve %%
contr_var_MIMO = [Msmin14,Msmin15,Msmin17,Msmin19,Msmin25,Msmin27,Msmin29,Msmin47;% bounds on S and T
                 KSmin14,KSmin15,KSmin17,KSmin19,KSmin25,KSmin27,KSmin29,KSmin47;% Bounds on KS
                 GammaSG14_min,GammaSG15_min,GammaSG17_min,GammaSG19_min,GammaSG25_min,GammaSG27_min,GammaSG29_min,GammaSG47_min;% Bounds on SG
                 KSGdmin141,KSGdmin151,KSGdmin171,KSGdmin191,KSGdmin251,KSGdmin271,KSGdmin291,KSGdmin471;% Bounds on KSGd1
                 KSGdmin142,KSGdmin152,KSGdmin172,KSGdmin192,KSGdmin252,KSGdmin272,KSGdmin292,KSGdmin472;% Bounds on KSGd2 
                 GammaSGd141_min,GammaSGd151_min,GammaSGd171_min,GammaSGd191_min,GammaSGd251_min,GammaSGd271_min,GammaSGd291_min,GammaSGd471_min;% Bounds on SGd1
                 GammaSGd142_min,GammaSGd152_min,GammaSGd172_min,GammaSGd192_min,GammaSGd252_min,GammaSGd272_min,GammaSGd292_min,GammaSGd472_min% Bounds on SGd2
    ]';

xlswrite('..\Hinf_Save\Controllability tables\Controllability_data_MIMO.xls',contr_var_MIMO,'Controllability_data_MIMO_1010','B2')
    
%%