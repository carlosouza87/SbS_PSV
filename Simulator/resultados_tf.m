close all
global iwaves1st iwavesmd iwind icontrol_psv ihdsuction imemory icoupling


%% Figure 1 - XY positions
% figure(1)
% lt = length(tsim);
% plot(variable.ship(1,1).eta(2,:),variable.ship(1,1).eta(1,:),'--')
% hold on
% L = data.ship(1,1).Lpp;
% B = data.ship(1,1).B;
% for kp1 = 1:round(lt/3):lt
%     [l1,l2,l3,l4,l5] = shipdraw(variable.ship(1,1).eta(1,kp1),variable.ship(1,1).eta(2,kp1),variable.ship(1,1).eta(6,kp1),L,B,1,-1,[0 0 1]);
% end
% plot(variable.ship(1,2).eta(2,:),variable.ship(1,2).eta(1,:),'r--')
% axis equal
% L = data.ship(1,2).Lpp;
% B = data.ship(1,2).B;
% for kp1 = 1:round(lt/3):lt
%     [l1,l2,l3,l4,l5] = shipdraw(variable.ship(1,2).eta(1,kp1),variable.ship(1,2).eta(2,kp1),variable.ship(1,2).eta(6,kp1),L,B,1,-1,[1 0 0]);
% end
% 
% clip
% 
% title('Posição do navios')
% xlabel('Posição em Y (m)')
% ylabel('Posição em X (m)')
% axis('equal')
% grid off
% 
% hold off


%% Figure 2 - CGs distance
figure(2)

for k1 = 1:lt
    d(k1) = abs(variable.ship(1,2).eta(2,k1) - variable.ship(1,1).eta(2,k1) - 0.5*(data.ship(1,1).B + data.ship(1,2).B));
end

plot(tsim,d)
title('Distância entre os costados','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Distância (m)','fontsize',12)
grid 
% 
%% Figure 3 - Hydrodynamic interaction loads

if ihdsuction == 1
    
figure(3)
subplot(2,1,1)
plot(tsim(1:lt-1),variable.ship(1).tau_hdsuction(2,1:lt-1))
title('Força de sucção em Y - FPSO','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Força (N)','fontsize',12)
grid minor

subplot(2,1,2)
plot(tsim(1:lt-1),variable.ship(2).tau_hdsuction(2,1:lt-1))
title('Força de sucção em Y - PSV','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Força (N)','fontsize',12)
grid minor

 end
%% Figure 4 - Movimento - FPSO
figure(4)
subplot(3,1,1)
plot(tsim,variable.ship(1).eta(1,1:lt))
title('Posição em X - FPSO','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('X (m)','fontsize',12)
grid minor

subplot(3,1,2)
plot(tsim,variable.ship(1).eta(2,1:lt))
title('Posição em Y - FPSO','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Y (m)','fontsize',12)
grid minor

subplot(3,1,3)
plot(tsim,variable.ship(1).eta(6,1:lt)*180/pi)
title('Posição em \Psi - FPSO','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('\Psi (graus)','fontsize',12)
grid minor
% 
% plot(tsim,variable.ship(1).eta(3,1:lt))
% title('Posição em Z - FPSO','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Z (m)','fontsize',12)
% grid minor

%% Figure 5 - Angular motions - ship 1
% figure(5)
% subplot(3,1,1)
% plot(tsim,variable.ship(1).eta(4,1:lt)*180/pi)
% title('Posição em \Phi - FPSO','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('\Phi (graus)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim,variable.ship(1).eta(5,1:lt)*180/pi)
% title('Posição em \Theta - FPSO','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('\Theta (graus)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim,variable.ship(1).eta(6,1:lt)*180/pi)
% title('Posição em \Psi - FPSO','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('\Psi (graus)','fontsize',12)
% grid minor

%% Figure 6 - Linear motions - ship 2
figure(6)
subplot(3,1,1)
plot(tsim,variable.ship(2).eta(1,1:lt))
title('Posição em X - PSV','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('X (m)','fontsize',12)
grid minor

subplot(3,1,2)
plot(tsim,variable.ship(2).eta(2,1:lt))
title('Posição em Y - PSV','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Y (m)','fontsize',12)
grid minor

subplot(3,1,3)
plot(tsim,variable.ship(2).eta(6,1:lt)*180/pi)
title('Posição em \Psi - PSV','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('\Psi (graus)','fontsize',12)
grid minor

% plot(tsim,variable.ship(2).eta(3,1:lt))
% title('Posição em Z - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Z (m)','fontsize',12)
% grid minor

%% Figure 7 - Angular motions - ship 2
% figure(7)
% subplot(3,1,1)
% plot(tsim,variable.ship(2).eta(4,1:lt)*180/pi)
% title('Posição em \Phi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('\Phi (graus)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim,variable.ship(2).eta(5,1:lt)*180/pi)
% title('Posição em \Theta - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('\Theta (graus)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim,variable.ship(2).eta(6,1:lt)*180/pi)
% title('Posição em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('\Psi (graus)','fontsize',12)
% grid minor

%% Figure 8 - Memory effects - ship 1 (linear motions)
% if imemory == 1
% figure(11)
% subplot(3,1,1)
% plot(tsim,variable.mu(1,1:lt))
% title('\mu - ship 1, x')
% xlabel('Time (s)')
% ylabel('\mu (N)')
% grid 
% subplot(3,1,2)
% plot(tsim,variable.mu(2,1:lt))
% title('\mu - ship 1, y')
% xlabel('Time (s)')
% ylabel('\mu (N)')
% grid 
% subplot(3,1,3)
% plot(tsim,variable.mu(3,1:lt))
% title('\mu - ship 1, z')
% xlabel('Time (s)')
% ylabel('\mu (N)')
% grid 
% end
%% Figure 9 - Memory effects - ship 2 (angular motions) %acho que é ship
%% 1.
% if imemory == 1
% figure(12)
% subplot(3,1,1)
% plot(tsim,variable.mu(4,1:lt))
% % title('\mu - ship 4, \phi') %essa linha estava originalmente e foi
% % substituída pela abaixo.
% title('\mu - ship 1, \phi')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% subplot(3,1,2)
% plot(tsim,variable.mu(5,1:lt))
% title('\mu - ship 1, \theta')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% subplot(3,1,3)
% plot(tsim,variable.mu(6,1:lt))
% title('\mu - ship 1, \psi')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% end
%% Figure 10 - Memory effects - ship 2 (linear motions)
% if imemory==1
% figure(13)
% subplot(3,1,1)
% plot(tsim,variable.mu(1,1:lt)) %VARIABLE.MU
% title('\mu - ship 2, x')
% xlabel('Time (s)')
% ylabel('\mu (N)')
% grid 
% subplot(3,1,2)
% plot(tsim,variable.mu(2,1:lt))
% title('\mu - ship 2, y')
% xlabel('Time (s)')
% ylabel('\mu (N)')
% grid 
% subplot(3,1,3)
% plot(tsim,variable.mu(3,1:lt))
% title('\mu - ship 2, z')
% xlabel('Time (s)')
% ylabel('\mu (N)')
% grid 
% end
%% Figure 11 - Memory effects - ship 2 (angular motions)
% if imemory==1
% figure(14)
% subplot(3,1,1)
% plot(tsim,variable.mu(4,1:lt))
% title('\mu - ship 2, \phi')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% subplot(3,1,2)
% plot(tsim,variable.mu(5,1:lt))
% title('\mu - ship 2, \theta')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% subplot(3,1,3)
% plot(tsim,variable.mu(6,1:lt))
% title('\mu - ship 2, \psi')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% % end
% %% Figure 15 - Memory effects - ship 2 (angular motions)- Matrix fixfreq
% figure(15)
% subplot(3,1,1)
% plot(tsim,variable.mu_fixfreq(4,1:lt))
% title('\mu_fixfreq - ship 2, \phi')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% subplot(3,1,2)
% plot(tsim,variable.mu_fixfreq(5,1:lt))
% title('\mu_fixfreq - ship 2, \theta')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 
% subplot(3,1,3)
% plot(tsim,variable.mu_fixfreq(6,1:lt))
% title('\mu_fixfreq - ship 2, \psi')
% xlabel('Time (s)')
% ylabel('\mu (Nm)')
% grid 

%% Figure 12 - Amarras loads 

lt1=lt-1;
tsim1=tsim(1:length(tsim)-1);

figure(12)
subplot(3,1,1)
plot(tsim1,variable.tau_amarras(1,1:lt1))
title('Tração nas amarras em X - FPSO','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Força (N)','fontsize',12)
grid minor
subplot(3,1,2)
plot(tsim1,variable.tau_amarras(2,1:lt1))
title('Tração das amarras em Y - FPSO','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Força (N)','fontsize',12)
grid minor
subplot(3,1,3)
plot(tsim1,variable.tau_amarras(6,1:lt1))
title('Momento das amarras em \Psi - FPSO','fontsize',13)
xlabel('Tempo (s)','fontsize',12)
ylabel('Momento (Nm)','fontsize',12)
grid minor

%% Figure 13 - forças externas e de controle 

% figure(13)
% subplot(3,1,1)
% plot(tsim1,variable.tau_ctr(7,1:lt1),tsim1,variable.tau_ext(7,1:lt1))
% title('Lei de controle e força externa em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid 
% 
% subplot(3,1,2)
% plot(tsim1,variable.tau_ctr(8,1:lt1),tsim1,variable.tau_ext(8,1:lt1))
% title('Lei de controle e força externa em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid 
% 
% subplot(3,1,3)
% plot(tsim1,variable.tau_ctr(12,1:lt1),tsim1,variable.tau_ext(12,1:lt1))
% title('Lei de controle e momento externo em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Momento (Nm)','fontsize',12)
% grid 

%% Figure 14 - Força estimada

% figure(14)
% subplot(3,1,1)
% plot(tsim1,variable.tau_b_hat_psv(7,1:lt1))
% title('Força estimada em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim1,variable.tau_b_hat_psv(8,1:lt1))
% title('Forca estimada em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim1,variable.tau_b_hat_psv(12,1:lt1))
% title('Momento estimado em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Momento (Nm)','fontsize',12)
% grid minor

%% Figure 15 - Força externa

% figure(15)
% subplot(3,1,1)
% plot(tsim1,variable.tau_ext(7,1:lt1))
% title('Força externa em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim1,variable.tau_ext(8,1:lt1))
% title('Forca externa em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim1,variable.tau_ext(12,1:lt1))
% title('Momento externo em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Momento (Nm)','fontsize',12)
% grid minor

%% Figure 16 - Posição estimada e real

% figure(16)
% subplot(3,1,1)
% plot(tsim1,variable.eta_hat(7,1:lt1),tsim1,variable.ship(2).eta(1,1:lt1))
% title('Posição estimada e real em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('X (m)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim1,variable.eta_hat(8,1:lt1),tsim1,variable.ship(2).eta(2,1:lt1))
% title('Posição estimada e real em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Y (m)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim1,variable.eta_hat(12,1:lt1)*180/pi,tsim1,variable.ship(2).eta(6,1:lt1)*180/pi)
% title('Aproamento estimado e real em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('\Psi (graus)','fontsize',12)
% grid minor

%% Figure 17 - Força de onda de primeira ordem

% if iwaves1st == 1
% figure(17)
% subplot(3,1,1)
% plot(tsim1,variable.ship(2).waves1st(1,1:lt1))
% title('Força de onda de 1º ordem em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim1,variable.ship(2).waves1st(2,1:lt1))
% title('Força de onda de 1º ordem em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim1,variable.ship(2).waves1st(6,1:lt1))
% title('Momento de onda de 1º ordem em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Momento (Nm)','fontsize',12)
% grid minor
% 
%  end
%% controle

%  if icontrol_psv == 1
%     
% figure(18)
% %plot(tsim1,variable.tau_ext(7,1:lt1),tsim1,variable.tau_b_hat_psv(7,1:lt1),tsim1, variable.tau_cont_aux(1,1:lt1))
% plot(tsim1,variable.tau_cont_aux(1,1:lt1), tsim1,variable.tau_cont_aux_Kp(1,1:lt1),tsim1, variable.tau_cont_aux_Kd(1,1:lt1))
% title('tau\_cont\_aux, tau\_cont\_aux\_Kp, tau\_cont\_aux\_Kd em X')
% grid
% 
% figure(19)
% %plot(tsim1,variable.tau_ext(8,1:lt1),tsim1,variable.tau_b_hat_psv(9,1:lt1), tsim1, variable.tau_cont_aux(2,1:lt1))
% plot(tsim1,variable.tau_cont_aux(2,1:lt1), tsim1,variable.tau_cont_aux_Kp(2,1:lt1),tsim1, variable.tau_cont_aux_Kd(2,1:lt1))
% title('tau\_cont\_aux, tau\_cont\_aux\_Kp, tau\_cont\_aux\_Kd em Y')
% grid
% 
% figure(20)
% %plot(tsim1,variable.tau_ext(12,1:lt1),tsim1,variable.tau_b_hat_psv(12,1:lt1),tsim1,variable.tau_cont_aux(3,1:lt1))
% plot(tsim1,variable.tau_cont_aux(3,1:lt1), tsim1,variable.tau_cont_aux_Kp(3,1:lt1),tsim1, variable.tau_cont_aux_Kd(3,1:lt1))
% title('tau\_cont\_aux - tau\_cont\_aux\_Kp - tau\_cont\_aux\_Kd em Z')
% grid
% 
%  end
%% Figure 22 - Força de deriva média

%  if iwavesmd == 1
% 
% figure(22)
% subplot(3,1,1)
% plot(tsim1,variable.ship(2).tau_wavesmd(7,1:lt1))
% title('Força de onda de deriva média em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim1,variable.ship(2).tau_wavesmd(8,1:lt1))
% title('Força de onda de deriva média em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim1,variable.ship(2).tau_wavesmd(12,1:lt1))
% title('Momento de onda de deriva média em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Momento (Nm)','fontsize',12)
% grid minor
% 
%  end

%%
% figure(23)
% plot(tsim1,variable.nu_hat(7,1:lt1),tsim1,variable.ship(2).nu(1,1:lt1))
% title('Velocidade real e estimada em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Velocidade (m/s)','fontsize',12)
% grid 
% 
% figure(24)
% plot(tsim1,variable.nu_hat(8,1:lt1),tsim1,variable.ship(2).nu(2,1:lt1))
% title('Velocidade real e estimada em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Velocidade (m/s)','fontsize',12)
% grid 
% 
% figure(25)
% plot(tsim1,variable.nu_hat(12,1:lt1),tsim1,variable.ship(2).nu(6,1:lt1))
% title('Velocidade real e estimada em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Velocidade (rad/s)','fontsize',12)
% grid 

%% mu
if imemory ==1

variable.t_plot = 0;
variable.t_plot(1:lt) = tsim(1:lt);


figure(26)
subplot(2,1,1)
plot(variable.t_plot,variable.mu(1,:)','-','LineWidth',2)
title('Função de memória em X - FPSO','fontsize',13)
xlabel ('Tempo (s)','fontsize',12)
ylabel ('Força (N)','fontsize',12)
grid on

subplot(2,1,2)
plot(variable.t_plot,variable.mu(2,:)','-','LineWidth',2)
title('Função de memória em Y - FPSO','fontsize',13)
xlabel ('Tempo (s)','fontsize',12)
ylabel ('Força (N)','fontsize',12)
grid on


figure(27)
subplot(2,1,1)
plot(variable.t_plot,variable.mu(7,:)','-','LineWidth',2)
title('Função de memória em X - PSV','fontsize',13)
xlabel ('Tempo (s)','fontsize',12)
ylabel ('Força (N)','fontsize',12)
grid on

subplot(2,1,2)
plot(variable.t_plot,variable.mu(8,:)','-','LineWidth',2)
title('Função de memória em Y - PSV','fontsize',13)
xlabel ('Tempo (s)','fontsize',12)
ylabel ('Força (N)','fontsize',12)
grid on

end

%% Figure 28 - Força da corrente

%  if icurr == 1
%     
% figure(28)
% subplot(3,1,1)
% plot(tsim(1:lt-1),variable.ship(2).tau_curr(7,:))
% title('Força de corrente em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim(1:lt-1),variable.ship(2).tau_curr(8,:))
% title('Força de corrente em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim(1:lt-1),variable.ship(2).tau_curr(12,:))
% title('Momento da corrente em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Momento (Nm)','fontsize',12)
% grid minor
% 
%  end
%% Figure 29 - Força do vento

%  if iwind == 1
%     
% figure(29)
% subplot(3,1,1)
% plot(tsim,variable.ship(2).tau_wind(7,:))
% title('Força do vento em X - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,2)
% plot(tsim,variable.ship(2).tau_wind(8,:))
% title('Força do vento em Y - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Força (N)','fontsize',12)
% grid minor
% 
% subplot(3,1,3)
% plot(tsim,variable.ship(2).tau_wind(12,:))
% title('Momento do vento em \Psi - PSV','fontsize',13)
% xlabel('Tempo (s)','fontsize',12)
% ylabel('Momento (Nm)','fontsize',12)
% grid minor
% 
%  end
%% Salvar gráficos
opcao = 0;

if opcao == 1
    
    figuras = [1 2 3 4 6 26 27];
    
    for i=1:length(figuras)

            fig = figuras(i);
            num_figura = num2str(fig);
            local = ['E:\TF\resultados\figura_' num_figura '.pdf'];
            position = [403 246 560 420];
            set(figure(fig),'Position',position);
            %-depsc color; -deps2 black and white; -depsc2 color 2; -dpdf pdf; -djpeg jpg
            print(figure(fig), '-dpdf', local);

    end

end

% 
% dist = [2 4 6 8 10 12 14 16 18 20 22 24];
% 
% forca = [1735 562 307 204 151 119 98 83 72 63 56 51];
% 
% 
% plot(dist,forca)
% title('Força de sucção em Y - PSV','fontsize',13)
% xlabel('Distância (m)','fontsize',12)
% ylabel('Força (kN)','fontsize',12)
% xlim([4 24])
% grid