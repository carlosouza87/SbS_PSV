function [tauc,Fip] = surgespeedctr(u,Uref,Uref_p,m,d,Kp,Ki,Fi,dt)
% [tauc] = surgespeedctr(u,Uref,Uref_p,m,d,Kp,Ki,Fi)
% 
% Calculates surge control force in order to keep the ship with a reference
% Uref forward speed, through de use of state feedback linearization
% method. Inputs are:
% u - current surge velocity [m/s]
% Uref - reference surge velocity [m/s]
% Uref_p - first time derivative of Uref_p [m/s^2]
% m - ship displacement [kg]
% d - surge quadratic drag coeficient (e.g. ITTC'57) []
% Kp - controller proportional gain [kg/s]
% Ki - controller integral gain [kg/m]
% Fi - accumulated integral parcel
% dt - time step [s]
% 
% Output is:
% tauc - surge control force [N]
% 
% Carlos Souza, 15/12/2010 - Universidade de São Paulo

tauc = m*(Uref_p-Kp*(u-Uref)+Fi) + d*abs(u)*u;
Fip = Ki*(u-Uref);
