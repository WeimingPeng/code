function sys = v4_new_4d_model1(x_in)

% Simplified model for riser slugging
% By: Esmaeil Jahanshahi
% August 2009, NTNU, Norway

% x1: Mass of gas in the pipelien (m_G1)
% x2: Mass of liquid in the pipeline (m_L1)
% x3: Mass of gas in the riser (m_G2)
% x4: Mass of liquid in the riser (m_L2)
global par

x = x_in(1:4);
u = x_in(5:7);
[xdot,y]=v4_new_4d_model0(x,u,par);

sys = [xdot;y];

end