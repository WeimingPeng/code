function [sys] = v1_new_4d_Observer(t,x_hat,u,y_m,output,par,ep)

% Observer based on the 4-state pipeline-riser model
% By: Esmaeil Jahanshahi
% May 2012, NTNU, Norway


% x1_hat: Mass of gas in the pipeline (m_G1)
% x2_hat: Mass of liquid in the pipeline (m_L1)
% x3_hat: Mass of gas in the riser (m_G2)
% x4_hat: Mass of liquid in the riser (m_L2)


[zdot,xhat] = v1_new_4d_Observer0(x_hat,u,y_m,par,ep);

if isequal(output,'derivatives')
    sys =zdot;
    if ~isreal(sys)
        disp('Complex')
        sys=0*sys;
    end
elseif isequal(output,'measurements')
    sys = xhat;
end


end
