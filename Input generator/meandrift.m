clear all;close all;clc

load waveforces

Hs = 2.5;
Tp = 9;

nfreqs = length(freqs);
nincid = length(unique_incid);

D01_1 = zeros(nfreqs,nincid);
D02_1 = zeros(nfreqs,nincid);
D06_1 = zeros(nfreqs,nincid);
Ph1_1 = zeros(nfreqs,nincid);
Ph2_1 = zeros(nfreqs,nincid);
Ph6_1 = zeros(nfreqs,nincid);
D01_2 = zeros(nfreqs,nincid);
D02_2 = zeros(nfreqs,nincid);
D06_2 = zeros(nfreqs,nincid);
Ph1_2 = zeros(nfreqs,nincid);
Ph2_2 = zeros(nfreqs,nincid);
Ph6_2 = zeros(nfreqs,nincid);

for k1 = 1:nfreqs
    for k2 = 1:nincid
        D01_1(k1,k2) = drift_amp1(k1,k2,1)*cos(drift_pha1(k1,k2,1)*pi/180);  % Wave drift coefficient - vlcc (surge)
        D01_2(k1,k2) = drift_amp2(k1,k2,1)*cos(drift_pha2(k1,k2,1)*pi/180);  % Wave drift coefficient - Suezmax (surge)
        D02_1(k1,k2) = drift_amp1(k1,k2,2)*cos(drift_pha1(k1,k2,2)*pi/180);  % Wave drift coefficient - vlcc (sway)
        D02_2(k1,k2) = drift_amp2(k1,k2,2)*cos(drift_pha2(k1,k2,2)*pi/180);  % Wave drift coefficient - Suezmax (sway)
        D06_1(k1,k2) = drift_amp1(k1,k2,6)*cos(drift_pha1(k1,k2,6)*pi/180);  % Wave drift coefficient - vlcc (yaw)
        D06_2(k1,k2) = drift_amp2(k1,k2,6)*cos(drift_pha2(k1,k2,6)*pi/180);  % Wave drift coefficient - Suezmax (yaw)]
        %         D01_1(k1,k2) = drift_amp1(k1,k2,1);  % Wave drift coefficient - vlcc (surge)
        %         D01_2(k1,k2) = drift_amp2(k1,k2,1);  % Wave drift coefficient - Suezmax (surge)
        %         D02_1(k1,k2) = drift_amp1(k1,k2,2);  % Wave drift coefficient - vlcc (sway)
        %         D02_2(k1,k2) = drift_amp2(k1,k2,2);  % Wave drift coefficient - Suezmax (sway)
        %         D06_1(k1,k2) = drift_amp1(k1,k2,6);  % Wave drift coefficient - vlcc (yaw)
        %         D06_2(k1,k2) = drift_amp2(k1,k2,6);  % Wave drift coefficient - Suezmax (yaw)
    end
end

w = linspace(freqs(1),freqs(nfreqs),50);
dw = w(2) - w(1);

for k1 = 1:nincid
    D01ext_1(:,k1) = interp1(freqs,D01_1(:,k1),w);
    D01ext_2(:,k1) = interp1(freqs,D01_2(:,k1),w);
    D02ext_1(:,k1) = interp1(freqs,D02_1(:,k1),w);
    D02ext_2(:,k1) = interp1(freqs,D02_2(:,k1),w);
    D06ext_1(:,k1) = interp1(freqs,D06_1(:,k1),w);
    D06ext_2(:,k1) = interp1(freqs,D06_2(:,k1),w);
end

spec = 2;

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

for k1 = 1:length(unique_incid)
    Fm1_1(k1,1) = dw*trapz(Sw.*D01ext_1(:,k1));
    Fm1_2(k1,1) = dw*trapz(Sw.*D01ext_2(:,k1));
    Fm2_1(k1,1) = -dw*trapz(Sw.*D02ext_1(:,k1));
    Fm2_2(k1,1) = -dw*trapz(Sw.*D02ext_2(:,k1));
    Fm6_1(k1,1) = -dw*trapz(Sw.*D06ext_1(:,k1));
    Fm6_2(k1,1) = -dw*trapz(Sw.*D06ext_2(:,k1));
end