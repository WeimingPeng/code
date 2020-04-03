function [De,Dd,Du] = scaling(z1,z2)
%SCALING Summary of this function goes here
%   Detailed explanation goes here
De=(diag([1 1 1 1 1 1 1 2 1 50 1]));
Dd=diag([0.05*0.36 0.05*8.64]);
%Dd=(diag([0.04 1]));
Du=diag([min(z1,1-z1) min(z2,1-z2)]);
end