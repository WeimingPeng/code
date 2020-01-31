function [sys] = v4_new_4d_model(t,x,u,output,par)

% Simplified model for riser slugging
% By: Esmaeil Jahanshahi
% August 2009, NTNU, Norway

% x1: Mass of gas in the pipelien (m_G1)
% x2: Mass of liquid in the pipeline (m_L1)
% x3: Mass of gas in the riser (m_G2)
% x4: Mass of liquid in the riser (m_L2)

[xdot,y]=v4_new_4d_model0(x,u,par);

if isequal(output,'derivatives')
    sys = xdot;
    if ~isreal(sys)
        disp('Complex')
        error('Complex')
        sys=0*sys;
    end
elseif isequal(output,'measurements')
    sys = y(1:3);

end 

end