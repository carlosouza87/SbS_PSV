function [angn] = norm0360(ang)
% [angn] = norm02pi(ang)
% Normaliza �ngulo entre 0� e 360�. Dado o �ngulo n�o-normalizado ang, a
% fun��o faz a normaliza��o e devolve o �ngulo normalizado angn.
% 
% Carlos Souza, 29/10/2009 - Universidade de S�o Paulo

if abs(ang) > 360
    ang = rem(ang,360);
end
if ang < 0
    ang = rem(ang,360);
    ang = 360 + ang;
end
angn = ang;