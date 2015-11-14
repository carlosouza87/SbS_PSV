function [X,Y,N] = hydrint_vantorre01(ksi,rho,L,B,T,U1,U2)

polCx = [0.0001 0.0001 -0.0016 -0.0029 0.0177 0.0320 -0.1209 -0.2143 ...
    0.5228 0.9021 -1.4304 -2.3700 2.4142 3.7079 -2.4103 -3.0982 1.3661 1.0409 -0.3988];
polCy = [0.0005 -0.0008 -0.0056 0.0087 0.0408 -0.0573 -0.1922 0.2351 ...
    0.5797 -0.5976 -1.0828 0.8950 1.1561 -0.7046 -0.5912 0.2179 0.0854];
polCn = [0.0001 0.0001 -0.0024 -0.0017 0.0275 0.0181 -0.1953 -0.1195 ...
    0.8824 0.4952 -2.5303 -1.2744 4.4370 1.9413 -4.3685 -1.5657 1.9953 0.5017 -0.2304];

if abs(ksi) <= 1
    Cx = polyval(polCx,ksi);
    Cy = polyval(polCy,ksi);
    Cn = polyval(polCn,ksi);
else
    Cx = 0;
    Cy = 0;
    Cn = 0;
end

X = Cx*(1/2)*rho*B*T*U1*U2;
Y = Cy*(1/2)*rho*B*T*U1*U2; %Y = Cy*(1/2)*rho*L*T*U1*U2; %a principio, como constatado só com mrn (=0), suction e fnd. A força de sway um pouco alta. 
N = Cn*(1/2)*rho*B*L*T*U1*U2; %talvez o escoamento só acelere no comprimento proporcional ao da segunda embarcação.

% if(t>11)
% end

% N = Cn/100*(1/2)*rho*B*L*T*U1*U2;
% N = Cn*(1/2)*rho*B*88*T*U1*U2;
%N = Cy*(1/2)*rho*B*L*T*U1*U2; %Cy de proposito
