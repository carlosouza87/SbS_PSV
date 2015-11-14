clear all;close all;clc

rho = 1025

load hydrodata_cs5d3m_D_3E5_old
% load hydrodata_recente

A22 = data.hydro.A22;
B22 = data.hydro.B22;
freqs = data.hydro.freqs;

for k1 = 1:92
    a11_rodrigo(k1) = A22(1,1,k1);
    a22_rodrigo(k1) = A22(2,2,k1);
    a66_rodrigo(k1) = A22(6,6,k1);
    
    b11_rodrigo(k1) = B22(1,1,k1);
    b22_rodrigo(k1) = B22(2,2,k1);
    b66_rodrigo(k1) = B22(6,6,k1);
end


load cnpq_10m
% load hydrodata_recente

A22 = data.hydro.A22;
B22 = data.hydro.B22;

for k1 = 1:92
    a11_carlos(k1) = A22(1,1,k1);
    a22_carlos(k1) = A22(2,2,k1);
    a66_carlos(k1) = A22(6,6,k1);
    
    b11_carlos(k1) = B22(1,1,k1);
    b22_carlos(k1) = B22(2,2,k1);
    b66_carlos(k1) = B22(6,6,k1);
end

figure(1)
subplot(3,2,1)
plot(freqs,a11_rodrigo)
hold on
plot(freqs,a11_carlos,'r')
grid on

subplot(3,2,2)
plot(freqs,b11_rodrigo)
hold on
plot(freqs,b11_carlos,'r')
grid on

subplot(3,2,3)
plot(freqs,a22_rodrigo)
hold on
plot(freqs,a22_carlos,'r')
grid on

subplot(3,2,4)
plot(freqs,b22_rodrigo)
hold on
plot(freqs,b22_carlos,'r')
grid on

subplot(3,2,5)
plot(freqs,a66_rodrigo)
hold on
plot(freqs,a66_carlos,'r')
grid on

subplot(3,2,6)
plot(freqs,b66_rodrigo)
hold on
plot(freqs,b66_carlos,'r')
grid on