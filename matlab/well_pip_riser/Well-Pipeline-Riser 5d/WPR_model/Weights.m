function [Wp1,Wp2,Wu,Wt] = Weights()
%WEIGHTS Summary of this function goes here
%   Detailed explanation goes here

% Weight without integral action
M1=2; wb1=0.05; a1=0.01;
Wp1 = 1/M1; %tf([1/M1 wb1],[1 wb1*a1]);

Wu = 10;

Wt1 = inv(tf([1/10 1],[1/0.1 1]));
Wt=Wt1*Wt1;

% Weight with integral action
M2=2; wb2=0.05; a2=0.01;
Wp2 =  tf([1/M2 wb2],[1 wb2*a2]);

% % Weight without integral action
% M1=2; wb1=0.01; a1=0.01;
% Wp1 = 1/M1; %tf([1/M1 wb1],[1 wb1*a1]);
% 
% % Weight with integral action
% M2=2; wb2=0.01; a2=0.01;
% Wp2 =  tf([1/M2 wb2],[1 wb2*a2]);
% 
% 
% Wu = 10;
% Wt1 = inv(tf([1/100 1],[1/0.2 1]));
% Wt=Wt1*Wt1;

end
