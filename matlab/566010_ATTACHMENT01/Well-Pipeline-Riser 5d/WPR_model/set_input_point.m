function u0 = set_input_point(a,b,c,d)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

u0(1) = a;
u0(2) = b;
u0(3) = c;
u0(4) = d;
save('input_point.mat','u0')

end

