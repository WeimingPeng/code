function [sys] = v2_new_5d_model(t,x,u,output,par)

% Simplified model for riser slugging
% By: Esmaeil Jahanshahi
% May 2011, NTNU, Norway

% x1: Total mass in the well (m_w)
% x2: Mass of gas in the pipeline (m_Gp)
% x3: Mass of liquid in the pipeline (m_Lp)
% x4: Mass of gas in the riser (m_Gr)
% x5: Mass of liquid in the riser (m_Lr)

[xdot,y]=v2_new_5d_model0(x,u,par);

if isequal(output,'derivatives')
    sys =xdot;
    if ~isreal(sys)
        disp('Complex')
        sys=0*sys;
    end
elseif isequal(output,'measurements')
    sys = [y(1)/1e5   %Y1 P_bh [bar] 
           y(2)/1e5   %Y2 P_wh [bar]
           y(3)       %Y3 w_mix_in  [kg/s]
           y(4)/1e5-1   %Y4 P_in [bar]
           y(5)/1e5-1   %Y5 P_rb [bar]
           y(6)/1e5   %Y6 DP_r [bar]
           y(7)/1e5   %Y7 P_rt [bar]
           y(8)       %Y8 Q_out [Litre/s]
           y(9)       %Y9 W_out [kg/s]
           y(10)      %Y10 Rho_t [kg/m3]
           y(11)      %Y11 Alpha_L [-]
           y(12)+0.5      %Y12 w_Gin [kg/s]
           y(13)+0.5];    %Y13 w_Lin [kg/s]
end
end
