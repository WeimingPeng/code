function [A,B,C,D] = v2_new_5d_linmodl_num(x0,u0,par)


A = zeros(5,5);
B = zeros(5,4);
C = zeros(11,5);
D = zeros(11,4);

dx = [0;0;0;0;0];
for i =1:5
    dx(i) = 1e-6;
    x = x0 + dx;
    xdot0 = v2_new_5d_model(0,x0,u0,'derivatives',par);
    xdot = v2_new_5d_model(0,x,u0,'derivatives',par);
    y0 = v2_new_5d_model(0,x0,u0,'measurements',par);
    y = v2_new_5d_model(0,x,u0,'measurements',par);
    A(:,i) = (xdot-xdot0)/dx(i);
    C(:,i) = (y(1:11)-y0(1:11))/dx(i);
    dx(i) = 0;
end

du = [0;0;0;0];
for i =1:4
    du(i) = 1e-6;
    u = u0 + du;
    xdot0 = v2_new_5d_model(0,x0,u0,'derivatives',par);
    xdot = v2_new_5d_model(0,x0,u,'derivatives',par);
    y0 = v2_new_5d_model(0,x0,u0,'measurements',par);
    y = v2_new_5d_model(0,x0,u,'measurements',par);
    B(:,i) = (xdot-xdot0)/du(i);
    D(:,i) = (y(1:11)-y0(1:11))/du(i);
    du(i) = 0;
end

end