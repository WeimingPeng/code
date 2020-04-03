clc
clear
x0 = 1;
y0 = 0;
r0 = [x0;y0];
tspan = [0, 15];
[t1,y1] = ode15s(@gfunc,tspan,r0);
figure(1)
plot(t1, y1(:,1));
hold on
plot(t1,y1(:,2));
hold off


function out = gfunc(t,y)
r1 = y(2)+t;
r2 = t;

out = [r1
    r2];
end