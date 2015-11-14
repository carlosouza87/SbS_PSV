function [A,K] = B2A(B,freqs,t,Ainf)

dt = t(2) - t(1);
tf = t(length(t));

w = linspace(freqs(1),freqs(length(freqs)),100);
dw = w(2) - w(1);
wf = w(length(w));

B = interp1(freqs,B,w);

for k1 = 1:length(t)
    for k2 = 1:length(w)
        integ(k2) = B(k2)*cos(w(k2)*t(k1));
    end
    K(k1) = 2/pi*dw*trapz(integ);
end
clear integ

for k1 = 1:length(w)
    for k2 = 1:length(t)
        integ(k2) = K(k2)*sin(w(k1)*t(k2));
    end
    A(k1,1) = Ainf - (1/(w(k1)+1e-10))*dt*trapz(integ);
end


% A = interp1(w,A,freqs);
% B = interp1(w,B,freqs);
% 
% for k1 = 1:length(freqs)
%     K(k1) = B(k1) + 1j*freqs(k1)*(A(k1)-Ainf);
% end

clear K

for k1 = 1:length(w)
    K(k1) = B(k1) + 1j*w(k1)*(A(k1)-Ainf);
end

A = interp1(w,A,freqs);

