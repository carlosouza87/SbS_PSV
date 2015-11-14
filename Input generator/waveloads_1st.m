function [Fwf1,Fwf2,Fwf3,Fwf4,Fwf5,Fwf6,Sl1,Sl2,Sl3,Sl4,Sl5,Sl6] = waveloads_1st(t,Hs,Tp,spec,amp,pha,freqs,randu)

lt = length(t);
nfreqs = length(freqs); 
nincid= size(amp,2);

w = linspace(freqs(1),freqs(nfreqs),nfreqs); 
dw = w(2) - w(1);

for k1 = 1:nincid
    amp1(:,k1) = interp1(freqs,amp(:,k1,1),w,'spline','extrap'); %checar se ela está dimensionalizada por rho*g no input generator
    amp2(:,k1) = interp1(freqs,amp(:,k1,2),w,'spline','extrap');
    amp3(:,k1) = interp1(freqs,amp(:,k1,3),w,'spline','extrap');
    amp4(:,k1) = interp1(freqs,amp(:,k1,4),w,'spline','extrap');
    amp5(:,k1) = interp1(freqs,amp(:,k1,5),w,'spline','extrap');
    amp6(:,k1) = interp1(freqs,amp(:,k1,6),w,'spline','extrap');
    pha1(:,k1) = interp1(freqs,pha(:,k1,1),w,'spline','extrap');
    pha2(:,k1) = interp1(freqs,pha(:,k1,2),w,'spline','extrap');
    pha3(:,k1) = interp1(freqs,pha(:,k1,3),w,'spline','extrap');
    pha4(:,k1) = interp1(freqs,pha(:,k1,4),w,'spline','extrap');
    pha5(:,k1) = interp1(freqs,pha(:,k1,5),w,'spline','extrap');
    pha6(:,k1) = interp1(freqs,pha(:,k1,6),w,'spline','extrap');
end

%% Wave spectrum
if spec == 1
    % Pierson-Moskowitz spectrum
    A = 4*pi^3*Hs^2/(0.710*Tp)^4;
    B = 16*pi^3/(0.710*Tp)^4;
    for k1 = 1:length(w)
        Sw(k1,1) = A*w(k1)^(-5)*exp(-B*w(k1)^(-4));
    end
elseif spec == 2
%     JONSWAP spectrum
    gamma = 1.1;
    w0 = 2*pi/Tp;
    %w0 = 1/Tp;
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

%% time-domain realization of the surge, sway and yaw motions
Fwf1 = zeros(lt,nincid);
Fwf2 = zeros(lt,nincid);
Fwf3 = zeros(lt,nincid);
Fwf4 = zeros(lt,nincid);
Fwf5 = zeros(lt,nincid);
Fwf6 = zeros(lt,nincid);

Sl1 = zeros(nfreqs,nincid);
Sl2 = zeros(nfreqs,nincid);
Sl3 = zeros(nfreqs,nincid);
Sl4 = zeros(nfreqs,nincid);
Sl5 = zeros(nfreqs,nincid);
Sl6 = zeros(nfreqs,nincid);

for k1 = 1:nincid
    
    for k2 = 1:nfreqs        %troquei 1:nfreqs por 1:50 %retroquei para o primeiro
        Sl1(k2,k1) = amp1(k2,k1)^2*Sw(k2); %função de transferencia ao quadrado vezes espectro da onda, dá o espectro da embarcação pela frequência
        Sl2(k2,k1) = amp2(k2,k1)^2*Sw(k2);
        Sl3(k2,k1) = amp3(k2,k1)^2*Sw(k2);
        Sl4(k2,k1) = amp4(k2,k1)^2*Sw(k2);
        Sl5(k2,k1) = amp5(k2,k1)^2*Sw(k2);
        Sl6(k2,k1) = amp6(k2,k1)^2*Sw(k2);
        A1 = sqrt(2*Sl1(k2,k1)*dw); %calculo de amplitude
        A2 = sqrt(2*Sl2(k2,k1)*dw);
        A3 = sqrt(2*Sl3(k2,k1)*dw);
        A4 = sqrt(2*Sl4(k2,k1)*dw);
        A5 = sqrt(2*Sl5(k2,k1)*dw);
        A6 = sqrt(2*Sl6(k2,k1)*dw);
        Fwf1(:,k1) = Fwf1(:,k1) + A1*cos(freqs(k2)*t'+pha1(k2,k1)+2*pi*randu(k2)); %aqui como uma somatória da contribuição de todas as frequencias
        Fwf2(:,k1) = Fwf2(:,k1) + A2*cos(freqs(k2)*t'+pha2(k2,k1)+2*pi*randu(k2));
        Fwf3(:,k1) = Fwf2(:,k1) + A3*cos(freqs(k2)*t'+pha3(k2,k1)+2*pi*randu(k2));
        Fwf4(:,k1) = Fwf2(:,k1) + A4*cos(freqs(k2)*t'+pha4(k2,k1)+2*pi*randu(k2));
        Fwf5(:,k1) = Fwf2(:,k1) + A5*cos(freqs(k2)*t'+pha5(k2,k1)+2*pi*randu(k2));
        Fwf6(:,k1) = Fwf6(:,k1) + A6*cos(freqs(k2)*t'+pha6(k2,k1)+2*pi*randu(k2));
  
        
% % % tau_waves1st
% % % tau_wavesmd
% % % tau_wind 
% % % tau_curr
% % % tau_ctr
% % % tau_hdsuction
% %  tau_fnd
% %  tau_mrn
% t
% u
% 
% % 
    end
end
% % faz a primeira vez pro FPSO e a segunda pro PSV
% % plot(freqs,Sl3,'o-','LineWidth',2)
% plot(t,Fwf3(:,1),'o-','LineWidth',2)
% xlabel('frequencies rad/s')hmmorish@usp.br	
% ylabel('Fwf3')
% % ylabel('Sw')
% grid on
Fwf2 = Fwf2;
Fwf4 = Fwf4;
Fwf6 = Fwf6;
