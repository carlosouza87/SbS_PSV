function [etap_hat,nup_hat,ksip_hat,bp_hat,eps] = nonl_observer(y,Jpsi,eta_hat,nu_hat,ksi_hat,b_hat,M,D,tauc,OMEGA,GAMMA,K1,K2,K3,K4,T,t)

%% Innovation calculation
y_hat = eta_hat + GAMMA*ksi_hat;
eps = y - y_hat;

%% Estimates derivatives
ksip_hat = OMEGA*ksi_hat + K1*eps;
etap_hat = Jpsi*nu_hat + K2*eps;
bp_hat = -T\b_hat + K3*eps;%eq. 3.21 do papper, T é uma matriz diagonal de constantes de tempo
%t;
nup_hat = M\(-D*nu_hat+Jpsi'*b_hat+tauc+Jpsi'*K4*eps);
%nup_hat = M\(-D*nu_hat+Jpsi'*b_hat+tauc+ K4*Jpsi'*eps); 

%% Reference
% ZAKARTCHOUK JUNIOR, A. - "Projeto de um observador passivo não-linear e de
% um controlador backstepping para navios de superfície". Universidade de
% São Paulo, 2010. Thesis (MSc)

