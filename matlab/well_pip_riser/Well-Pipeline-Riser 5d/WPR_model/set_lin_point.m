function [] = set_lin_point(a,b)
%UNTITLED 11 Summary of this function goes here
%   Detailed explanation goes here
output = 1;
z(1) = a;
z(2) = b;

save('lin_point.mat','z')

end

