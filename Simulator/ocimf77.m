function [Xw,Yw,Nw] = ocimf77(coef,gamma,rho_wind,Vrel,At,Al,L,bow,load)
% [Xw,Yw,Nw] = ocimf77(coef,gamma,rho_wind,Vrel,At,Al,L,bow,load)
% 
% Calculates forces and moments (horizontal DOF) due to wind through the 
% method proposed in (OCIMF, 1977). The inputs are:
% coef - matrix with the drag coefficients, like presented below
% gamma - the incidence angle of the wind in the hull [0º - 360º]
% rho_wind - air density [kg/m^3]
% Vrel - velocity of the air related to the ship [m/s]
% At - transverse projected area (above surface) [m^2]
% Al - lateral projected area (above surface) [m^2]
% L - length overall [m]
% bow - bow type (0 for conventional, 1 for bulbous)
% load - loading condition [0 for ballasted, 1 for loaded)
% 
% The outputs are:
%
% Xw - wind force, surge [N]
% Yw - wind force, sway [N]
% Nw - wind moment, yaw [Nm]
% 
% Each line of the "coef" matrix must have the following configuration
% 
% gamma Cxl Cxbbb Cxbcb Cyl Cyb Cnl Cnb
% 
% where:
% Cxl - coefficient for surge direction (loaded ship)
% Cxbbb - coefficient for surge direction (ballasted ship, bulbous bow)
% Cxblb - coefficient for surge direction (ballasted ship, conventional bow)
% Cyl - coefficient for sway direction (loaded ship)
% Cyb - coefficient for sway direction (ballasted ship)
% Cnl - coefficient for yaw direction (loaded ship)
% Cnb - coefficient for yaw direction (ballasted ship)
% 
% Carlos Souza, 20/12/2010 - Universidade de São Paulo

%% 'gamma', 'bow' and 'load' values testing
if gamma < 0 || gamma > 360
    error('Variable "gamma" must be in the range 0º - 360º.')
elseif gamma <= 180
    sinal = 1;
elseif gamma > 180
    gamma = 360 - gamma;
    sinal = -1;
end
% gamma
% sinal
ang = coef(:,1);
if load == 0
    if bow == 0
        Cx = interp1(ang,coef(:,4),gamma,'linear','extrap');
    elseif bow == 1
        Cx = interp1(ang,coef(:,3),gamma,'linear','extrap');
    else
        error('Variable "bow" must assume values 0 or 1.')
    end
    Cy = sinal*interp1(ang,coef(:,6),gamma,'linear','extrap');
    Cn = sinal*interp1(ang,coef(:,8),gamma,'linear','extrap');
    %disp('entrei')
elseif load == 1;
    Cx = sinal*interp1(ang,coef(:,2),gamma,'linear','extrap');
    Cy = sinal*interp1(ang,coef(:,5),gamma,'linear','extrap');
    Cn = sinal*interp1(ang,coef(:,7),gamma,'linear','extrap');
else
    error('Variable "load" must assume values 0 or 1.')
end

%% wind loads calculation
Xw = (1/7.6)*Cx*rho_wind*Vrel^2*At;
Yw = (1/7.6)*Cy*rho_wind*Vrel^2*Al;
Nw = (1/7.6)*Cn*rho_wind*Vrel^2*Al*L;

%% Reference
% OCIMF - "Prediction of Wind and Current Loads on VLCCs." London, U.K.: Oil
% Companies International Maritime Forum, 1977. 1-77.