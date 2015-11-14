function [Sw] = PMmod(w,Hs,Tp)
% [Sw] = PMmod(w,Hs,Tp)
% Calcula o espectro de Pierson-Moskowitz modificado de acordo com a ITTC
% (1978). Dado um vetor com as frequências discretizadas w [rad/s], a
% altura signifactiva de ondas Hs [m] e o período modal Tp [s], a função
% calcula o vetor Sw [m.s^2], com os valores do espectro em função de w.
% 
% Carlos Souza, 12/11/2009 - Universidade de Sâo Paulo

lw = length(w);
A = 4*pi^3*Hs^2/(0.710*Tp)^4;
B = 16*pi^3/(0.710*Tp)^4;
if w(1) == 0
    w(1) = w(2);
end

Sw = zeros(lw,1);
for k1 = 1:lw
    Sw(k1,1) = A*w(k1)^(-5)*exp(-B*w(k1)^(-4));
end