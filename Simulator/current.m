function [Xc,Yc,Nc] = current(Lpp,B,T,S,Cb,Cy,lL,ur,vr,r)
rho = 1025;
nu = 1e-6;

%% Reynolds number calculation
Ur = sqrt(ur^2+vr^2); % current velocity related to the ship hull
Re = Ur*Lpp/nu + 1000; % Reynolds number []
if Re == 0
    Re = 10;
end

%% coefficients calculation
if Ur == 0
    Xc = 0;
    Yc = 0;
    Nc = 0;
else
    alpha = atan2(-vr,-ur);
    C1 = (0.09375*S/((log10(Re)-2)^2*T*Lpp))*cos(alpha) + (pi*T/(8*Lpp))*(cos(3*alpha)-cos(alpha));
    C2 = (Cy-pi*T/(2*Lpp))*sin(alpha)*abs(sin(alpha)) + pi*T/(2*Lpp)*sin(alpha)^3 ...
        + pi*T/Lpp*(1+0.4*Cb*B/T)*sin(alpha)*abs(cos(alpha));
    C6 = -lL*(Cy-pi*T/(2*Lpp))*sin(alpha)*abs(sin(alpha)) - pi*T/Lpp*sin(alpha)*cos(alpha)...
        -((1+abs(cos(alpha)))/2)^2*pi*T/Lpp*(0.5-2.4*T/Lpp)*sin(alpha)*abs(cos(alpha));
    Cd2 = pi*T/(2*Lpp)*(1-4.4*B/Lpp+0.16*B/T);
    Cd6 = pi*T/(4*Lpp)*(1+0.16*B/T-2.2*B/Lpp);
    
    %% current loads
    Xc = 0.5*rho*T*Lpp*C1*Ur^2;
    Yc = 0.5*rho*T*Lpp*C2*Ur^2;
    Nc = 0.5*rho*T*Lpp^2*C6*Ur^2;
    
    %% damping due to yaw
    Xd = -1/4*rho*pi*T^2*Lpp*vr*r - 1/16*rho*pi*T^2*Lpp^2*sign(ur)*r^2;
    Yd = 1/2*rho*T*Lpp^2*Cd2*ur*r - 0.035*rho*T*Lpp^2*vr*r - 0.007*rho*T*Lpp^3*abs(r)*r;
    Nd = -1/2*rho*T*Lpp^3*Cd6*abs(ur)*r-3/20*rho*T*Lpp^3*Cy*abs(vr)*r-1/32*rho*T*Lpp^4*Cy*abs(r)*r;
    
    %% total current loads
    Xc = Xc + Xd; %Xc é a maior carga que influencia, calculada ora pro FPSO ora pro PSV.
    Yc = Yc + Yd;
    Nc = Nc + Nd;
end

%% Reference
% LEITE, A.J.P., ARANHA, J.A.P., UMEDA, C., de CONTI, M.B. (1998) - "Current
% forces in tankers and bifurcation of equilibrium of systems: hydrodynamic
% model and experiments". Applied Ocean Research, 20, 145-156, Elsevier
% Science