function [Fm1,Fm2,Fm6] = meandrift(Hs,Tp,spec,amp,pha,freqs) 

nfreqs = length(freqs);
nincid = size(amp,2);
% nincid
% amp

D01 = zeros(nfreqs,nincid);
D02 = zeros(nfreqs,nincid);
D06 = zeros(nfreqs,nincid);
D01ext = zeros(50,nincid);
D02ext = zeros(50,nincid);
D06ext = zeros(50,nincid);

for k1 = 1:nfreqs
    for k2 = 1:nincid
        D01(k1,k2) = amp(k1,k2,1)*cos(pha(k1,k2,1)*pi/180);  % Wave drift coefficient (surge)
        D02(k1,k2) = amp(k1,k2,2)*cos(pha(k1,k2,2)*pi/180);  % Wave drift coefficient (sway)        
        D06(k1,k2) = amp(k1,k2,6)*cos(pha(k1,k2,6)*pi/180);  % Wave drift coefficient (yaw)
    end
end

w = linspace(freqs(1),freqs(nfreqs),50);
dw = w(2) - w(1);

for k1 = 1:nincid
    D01ext(:,k1) = interp1(freqs,D01(:,k1),w,'spline','extrap');
    D02ext(:,k1) = interp1(freqs,D02(:,k1),w,'spline','extrap');
    D06ext(:,k1) = interp1(freqs,D06(:,k1),w,'spline','extrap');
end

Sw = zeros(length(w),1);
if spec == 1
    % Pierson-Moskowitz spectrum
    A = 4*pi^3*Hs^2/(0.710*Tp)^4;
    B = 16*pi^3/(0.710*Tp)^4;
    for k1 = 1:length(w)
        Sw(k1,1) = A*w(k1)^(-5)*exp(-B*w(k1)^(-4));
    end
elseif spec == 2
    % JONSWAP spectrum
    gamma = 1.1;
    w0 = 2*pi/Tp;
    for k1 = 1:length(w)
        if w(k1) <= w0
            sigma = 0.07;
        else
            sigma = 0.09;
        end
        Y = exp(-(w(k1)-w0)^2/(2*sigma^2*w0^2));        
        Sw(k1,1) = ((5*Hs^2*w0^4*(1-0.287*log(gamma)))/(16*w(k1)^5))*exp(-(5/4)*(w0/w(k1))^4)*gamma^Y;
    end
end
%nincid é 13
Fm1 = zeros(nincid,1);
Fm2 = zeros(nincid,1);
Fm6 = zeros(nincid,1);
for k1 = 1:nincid % linha 4 daqui ... nincid = size(amp,2);
    Fm1(k1,1) = 2*dw*trapz(Sw.*D01ext(:,k1));
    Fm2(k1,1) = 2*dw*trapz(Sw.*D02ext(:,k1));
    Fm6(k1,1) = 2*dw*trapz(Sw.*D06ext(:,k1));
end
% Fm1
% Fm2
% Fm6