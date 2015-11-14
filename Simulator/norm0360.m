function [angn] = norm0360(ang)
% [angn] = norm02pi(ang)
% Normaliza ângulo entre 0º e 360º. Dado o ângulo não-normalizado ang, a
% função faz a normalização e devolve o ângulo normalizado angn.
% 
% Carlos Souza, 29/10/2009 - Universidade de São Paulo

if abs(ang) > 360
    ang = rem(ang,360);
end
if ang < 0
    ang = rem(ang,360);
    ang = 360 + ang;
end
angn = ang;