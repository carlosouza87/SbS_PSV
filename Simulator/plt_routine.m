close all
global  iplot_control

%% Figure 1 - XY positions
figure(1)
lt = length(tsim);
plot(variable.ship(1,1).eta(2,:),variable.ship(1,1).eta(1,:),'--')
hold on
L = data.ship(1,1).Lpp;
B = data.ship(1,1).B;
for kp1 = 1:round(lt/3):lt
    [l1,l2,l3,l4,l5] = shipdraw(variable.ship(1,1).eta(1,kp1),variable.ship(1,1).eta(2,kp1),variable.ship(1,1).eta(6,kp1),L,B,1,-1,[0 0 1]);
end
plot(variable.ship(1,2).eta(2,:),variable.ship(1,2).eta(1,:),'r--')
axis equal
L = data.ship(1,2).Lpp;
B = data.ship(1,2).B;
for kp1 = 1:round(lt/3):lt
    [l1,l2,l3,l4,l5] = shipdraw(variable.ship(1,2).eta(1,kp1),variable.ship(1,2).eta(2,kp1),variable.ship(1,2).eta(6,kp1),L,B,1,-1,[1 0 0]);
end

%  clip

title('Ship position')
xlabel('Y axis (m)')
ylabel('X axis (m)')
axis('equal')
grid on

hold off

%% Figure 2 - CGs distance
figure(2)
title('CGs distance')
for k1 = 1:lt
    d(k1) = sqrt((variable.ship(1,1).eta(1,k1)-variable.ship(1,2).eta(1,k1))^2+(variable.ship(1,1).eta(2,k1)-variable.ship(1,2).eta(2,k1))^2);
end
plot(tsim,d)
xlabel('Time (s)')
ylabel('Distance (m)')
grid minor
% 
% Figure 3 - Hydrodynamic interaction loads
figure(3)
subplot(2,2,1)
plot(tsim(1:lt-1),variable.ship(1).tau_hdsuction(2,1:lt-1))
title('ship 1 - Hydrodynamic interaction sway force')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(2,2,2)
plot(tsim(1:lt-1),variable.ship(1).tau_hdsuction(6,1:lt-1))
title('ship 1 - Hydrodynamic interaction yaw moment')
xlabel('Time (s)')
ylabel('Moment (Nm)')
grid minor
subplot(2,2,3)
plot(tsim(1:lt-1),variable.ship(2).tau_hdsuction(2,1:lt-1))
title('ship 2 - Hydrodynamic interaction sway force')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(2,2,4)
plot(tsim(1:lt-1),variable.ship(2).tau_hdsuction(6,1:lt-1))
title('ship 2 - Hydrodynamic interaction yaw moment')
xlabel('Time (s)')
ylabel('Moment (Nm)')
grid minor


%% Figure 7 - Linear motions - ship 1
figure(7)
subplot(3,1,1)
plot(tsim,variable.ship(1).eta(1,1:lt))
title('X motions - ship 1')
xlabel('Time (s)')
ylabel('X (m)')
grid minor
subplot(3,1,2)
plot(tsim,variable.ship(1).eta(2,1:lt))
title('Y motions - ship 1')
xlabel('Time (s)')
ylabel('Y (m)')
grid minor
subplot(3,1,3)
plot(tsim,variable.ship(1).eta(3,1:lt))
title('Z motions - ship 1')
xlabel('Time (s)')
ylabel('Z (m)')
grid minor

%% Figure 8 - Angular motions - ship 1
figure(8)
subplot(3,1,1)
plot(tsim,variable.ship(1).eta(4,1:lt)*180/pi)
title('\Phi motions - ship 1')
xlabel('Time (s)')
ylabel('\Phi (deg)')
grid minor
subplot(3,1,2)
plot(tsim,variable.ship(1).eta(5,1:lt)*180/pi)
title('\Theta motions - ship 1')
xlabel('Time (s)')
ylabel('\Theta (deg)')
grid minor
subplot(3,1,3)
plot(tsim,variable.ship(1).eta(6,1:lt)*180/pi)
title('\Psi motions - ship 1')
xlabel('Time (s)')
ylabel('\Psi (deg)')
grid minor

%% Figure 9 - Linear motions - ship 2
figure(9)
subplot(3,1,1)
plot(tsim,variable.ship(2).eta(1,1:lt))
title('X motions - ship 2')
xlabel('Time (s)')
ylabel('X (m)')
grid minor
subplot(3,1,2)
plot(tsim,variable.ship(2).eta(2,1:lt))
title('Y motions - ship 2')
xlabel('Time (s)')
ylabel('Y (m)')
grid minor
subplot(3,1,3)
plot(tsim,variable.ship(2).eta(3,1:lt))
title('Z motions - ship 2')
xlabel('Time (s)')
ylabel('Z (m)')
grid minor

%% Figure 10 - Angular motions - ship 2
figure(10)
subplot(3,1,1)
plot(tsim,variable.ship(2).eta(4,1:lt)*180/pi)
title('\Phi motions - ship 2')
xlabel('Time (s)')
ylabel('\Phi (deg)')
grid minor
subplot(3,1,2)
plot(tsim,variable.ship(2).eta(5,1:lt)*180/pi)
title('\Theta motions - ship 2')
xlabel('Time (s)')
ylabel('\Theta (deg)')
grid minor
subplot(3,1,3)
plot(tsim,variable.ship(2).eta(6,1:lt)*180/pi)
title('\Psi motions - ship 2')
xlabel('Time (s)')
ylabel('\Psi (deg)')
grid minor


%% Figure 16 - Amarras loads 
lt1=lt-1;
tsim1=tsim(1:length(tsim)-1);
figure(16)
subplot(3,1,1)
plot(tsim1,variable.tau_amarras(1,1:lt1))
title('Tração Amarra em X')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,2)
plot(tsim1,variable.tau_amarras(2,1:lt1))
title('Tração Amarra em Y')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,3)
plot(tsim1,variable.tau_amarras(6,1:lt1))
title('Momento Amarra em torno de Z')
xlabel('Time (s)')
ylabel('Momentum (Nm)')
grid minor

%% Figure 17 - forças externas
% variable.tau_ctr(:,ktime) = tau_ctr;
%   variable.tau_b_hat_psv(:,ktime)=tau_b_hat;
if iplot_control == 1

    figure(17)
subplot(3,1,1)
plot(tsim1,variable.tau_ctr(7,1:lt1),tsim1,variable.tau_ext(7,1:lt1))
title('lei de controle e força externa em X')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,2)
plot(tsim1,variable.tau_ctr(8,1:lt1),tsim1,variable.tau_ext(8,1:lt1))
title('Lei de controle e força externa em Y')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,3)
plot(tsim1,variable.tau_ctr(12,1:lt1),tsim1,variable.tau_ext(12,1:lt1))
title('Lei de controle e momento externo em torno de Z')
xlabel('Time (s)')
ylabel('Momentum (Nm)')
grid minor


figure(18)
subplot(3,1,1)
plot(tsim1,variable.tau_b_hat_psv(7,1:lt1))
title('força estimada em X')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,2)
plot(tsim1,variable.tau_b_hat_psv(8,1:lt1))
title('forca estimada em Y')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,3)
plot(tsim1,variable.tau_b_hat_psv(12,1:lt1))
title('momento estimado em torno de Z')
xlabel('Time (s)')
ylabel('Momentum (Nm)')
grid minor

figure(19)
subplot(3,1,1)
plot(tsim1,variable.tau_ext(7,1:lt1))
title('força externa em X')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,2)
plot(tsim1,variable.tau_ext(8,1:lt1))
title('forca externa em Y')
xlabel('Time (s)')
ylabel('Force (N)')
grid minor
subplot(3,1,3)
plot(tsim1,variable.tau_ext(12,1:lt1))
title('momento externo em torno de Z')
xlabel('Time (s)')
ylabel('Momentum (Nm)')
grid minor


figure(20)
subplot(3,1,1)
plot(tsim1,variable.eta_hat(7,1:lt1),tsim1,variable.ship(2).eta(1,1:lt1))
title('X estimado e real')
xlabel('Time (s)')
ylabel('m)')
grid minor
subplot(3,1,2)
plot(tsim1,variable.eta_hat(8,1:lt1),tsim1,variable.ship(2).eta(2,1:lt1))
title('Y estimado e real')
xlabel('Time (s)')
ylabel('m')
grid minor
subplot(3,1,3)
plot(tsim1,variable.eta_hat(12,1:lt1)*180/pi,tsim1,variable.ship(2).eta(6,1:lt1)*180/pi)
title('angulo de aproamento e real')
xlabel('Time (s)')
ylabel('graus')
grid minor

figure(21)
subplot(3,1,1)
plot(tsim1,variable.ship(2).waves1st(1,1:lt1))
title('força primeira ordem X')
xlabel('Time (s)')
ylabel('N)')
grid minor
subplot(3,1,2)
plot(tsim1,variable.ship(2).waves1st(2,1:lt1))
title('força de primeira ordem Y')
xlabel('Time (s)')
ylabel('N')
grid minor
subplot(3,1,3)
plot(tsim1,variable.ship(2).waves1st(6,1:lt1)*180/pi)
title('momento de primeira ordem')
xlabel('Time (s)')
ylabel('Nm')
grid minor




figure(22)
%plot(tsim1,variable.tau_ext(7,1:lt1),tsim1,variable.tau_b_hat_psv(7,1:lt1),tsim1, variable.tau_cont_aux(1,1:lt1))
plot(tsim1,variable.tau_cont_aux(1,1:lt1), tsim1,variable.tau_cont_aux_Kp(1,1:lt1),tsim1, variable.tau_cont_aux_Kd(1,1:lt1))

title('tau_cont_aux, tau_cont_aux_Kp, tau_cont_aux_Kd em X')

figure(23)
%plot(tsim1,variable.tau_ext(8,1:lt1),tsim1,variable.tau_b_hat_psv(9,1:lt1), tsim1, variable.tau_cont_aux(2,1:lt1))
plot(tsim1,variable.tau_cont_aux(2,1:lt1), tsim1,variable.tau_cont_aux_Kp(2,1:lt1),tsim1, variable.tau_cont_aux_Kd(2,1:lt1))
title('tau_cont_aux, tau_cont_aux_Kp, tau_cont_aux_Kd em Y')

figure(24)
%plot(tsim1,variable.tau_ext(12,1:lt1),tsim1,variable.tau_b_hat_psv(12,1:lt1),tsim1,variable.tau_cont_aux(3,1:lt1))
plot(tsim1,variable.tau_cont_aux(3,1:lt1), tsim1,variable.tau_cont_aux_Kp(3,1:lt1),tsim1, variable.tau_cont_aux_Kd(3,1:lt1))
title('tau_cont_aux, tau_cont_aux_Kp, tau_cont_aux_Kd em torno de Z')


figure(25)
plot(tsim1,variable.tau_ext(7,1:lt1),tsim1,variable.ship(2).tau_curr(7,1:lt1),tsim1, variable.ship(2).tau_wavesmd(7,1:lt1))
title('força externa, corrente e wavesmod X')

figure(26)
plot(tsim1,variable.tau_ext(8,1:lt1),tsim1,variable.ship(2).tau_curr(8,1:lt1), tsim1, variable.ship(2).tau_wavesmd(8,1:lt1))
title('força externa e corrente e wavesmod em Y')

figure(27)
plot(tsim1,variable.tau_ext(12,1:lt1),tsim1,variable.ship(2).tau_curr(12,1:lt1),tsim1,variable.ship(2).tau_wavesmd(12,1:lt1))
title('momento extern0 e corrente e wavesmod em torno de Z')

%variable.ship(k1).tau_curr(:,ktime) = tau_curr;

figure(28)
plot(tsim1,variable.ship(2).nu(1,1:lt1),tsim1,variable.nu_hat(7,1:lt1))
title('velocidade real e estimada X')

figure(29)
plot(tsim1,variable.ship(2).nu(2,1:lt1),tsim1,variable.nu_hat(8,1:lt1))
title('velocidade real e estimada Y')

figure(30)
plot(tsim1,variable.ship(2).nu(6,1:lt1),tsim1,variable.nu_hat(12,1:lt1))
title('velocidade angular e estimada em torno de Z')

end

%variable.ship(2).nu(:,ktime) = y(19:24,ktime);


