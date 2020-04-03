%% function:  slug_mode0
% university: SWPU
% time: 2020
% design by WeimingPeng

%% main code
function [sys] =a_model_slug(t,x,u,output,par)

% x1: Mass of gas in the pipelien (m_G1)
% x2: Mass of liquid in the pipeline (m_L1)
% x3: Mass of gas in the riser (m_G2)
% x4: Mass of liquid in the riser (m_L2)

[xdot,y]=a_model0_slug(x,u,par);

if isequal(output,'derivatives')
    sys = xdot;
    if ~isreal(sys)
        disp('Complex')
        sys=0*sys;
    end
elseif isequal(output,'measurements')
    sys = y(1:3);
end 
end