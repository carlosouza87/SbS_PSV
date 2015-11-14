function [ang_norm] = norm02pi(ang)

if abs(ang) > 2*pi
    ang = rem(ang,2*pi);
end

if ang < 0
    ang = rem(ang,2*pi);
    ang = 2*pi + ang;
end
ang_norm = ang;