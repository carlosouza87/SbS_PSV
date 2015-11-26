%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%     Simulation of side-by-side operations between a PSV and an FPSO    %%

clear variables global;close all;clc

global ktime dt data variable
global iwaves1st iwavesmd iwind icontrol_psv ihdsuction imemory icoupling iplot_control ivalida icurr imooring
t1 = cputime;
nn=0;

%% Time parameters
dt = .1;% time step [s]   
tfinal = 100;% total simulation time [s]
lt = tfinal/dt+1;   % number of time steps including zero []
tsim = 0:dt:tfinal; % time vector
ktime = 1;          % time step [s]
data.constants.dt = dt;
variable.dt = dt;
variable.lt=lt;

%% Simulation type and number of DOFs
% Determine type of simulation, "isimtype": 
% isimtype = 1, for Cummins equation 
% isimtype = 2, for LF + WF superposition 
isimtype = 1;   

% Determine number o degrees of freedom, "idof":
% idof = 1, for 6 DOF
% idof = 2, for 3 DOF
idof = 1;

%% Flags for switching modules on/off
iwaves1st = 1;       % Flag for 1st order wave loads
iwavesmd = 1;        % Flag for meand drift loads
iwind = 0;           % Flag for wind loads
icurr = 0;           % Flag for current loads
icontrol_psv = 0;    % Flag for the psv control system
ihdsuction = 0;      % Flag for hydrodynamci suction loads
icoupling = 0;       % Flag for coupled dynamics
imooring = 0;        % Flag for FPSO mooring system
% imemory = 0;       % Flag for memory effects
% ivalida = 0;       % Flag to confirm validation of convolution evaluation

%% Loading of simulation data
dimensions      % m-file with definition of ships dimensions for each loading case
simdata         % Read simulation parameters and organize simulation data into structures

%% Environmental loads
% Waves
load randu1.txt -ascii      % list with uniformly distributed values
spec = 2;                   % flag for wave spectrum (1 = Pierson-Moskowitz, 2 = JONSWAP)
data.ship(1).waveincid = [data.ship(1).waveincid;360];
data.ship(2).waveincid = [data.ship(2).waveincid;360];

% 1st order wave loads calculation
if iwaves1st == 1
    for k1 = 1:2
        amp = data.ship(k1).FTF_amp;
        pha = data.ship(k1).FTF_pha;
        freqs = data.ship(k1).wavefreqs;
        [Fwf1,Fwf2,Fwf3,Fwf4,Fwf5,Fwf6,Sl1,Sl2,Sl3,Sl4,Sl5,Sl6] = waveloads_1st(tsim,Hs,Tp,spec,amp,pha,freqs,randu1);
        data.ship(k1).Fwf(:,:,1) = [Fwf1 Fwf1(:,1)];
        data.ship(k1).Fwf(:,:,2) = [Fwf2 Fwf2(:,1)];
        data.ship(k1).Fwf(:,:,3) = [Fwf3 Fwf3(:,1)];
        data.ship(k1).Fwf(:,:,4) = [Fwf4 Fwf4(:,1)];
        data.ship(k1).Fwf(:,:,5) = [Fwf5 Fwf5(:,1)];
        data.ship(k1).Fwf(:,:,6) = [Fwf6 Fwf6(:,1)];
        data.ship(k1).Sl(:,:,1) = [Sl1 Sl1(:,1)];
        data.ship(k1).Sl(:,:,2) = [Sl2 Sl2(:,1)];
        data.ship(k1).Sl(:,:,3) = [Sl3 Sl3(:,1)];
        data.ship(k1).Sl(:,:,4) = [Sl4 Sl4(:,1)];
        data.ship(k1).Sl(:,:,5) = [Sl5 Sl5(:,1)];
        data.ship(k1).Sl(:,:,6) = [Sl6 Sl6(:,1)];
    end
end

% 2nd order wave loads calculation
if iwavesmd == 1
    for k1 = 1:2
        amp = data.ship(k1).drift_amp; %(ship_matrices)data.ship(k1).drift_amp = drift_amp1 -> drift_amp1(nfreqs-k1+1,k2,k3) = amp(cont+k3,1)-> amp* veio do arquivo .9
        pha = data.ship(k1).drift_pha;
        freqs = data.ship(k1).wavefreqs;
        [Fm1,Fm2,Fm6] = meandrift(Hs,Tp,spec,amp,pha,freqs); % meandrift(Hs,Tp,spec,amp,pha,freqs) -> [Fm1,Fm2,Fm6]-> data.ship(k1).Fmd (de 3 colunas)
        data.ship(k1).Fmd(:,1) = [Fm1;Fm1(1)]; %recebe o vetor todo e cria uma ultima linha com o Fm1(1)
        data.ship(k1).Fmd(:,2) = [Fm2;Fm2(1)]; %quando olha-se a variavel pelo �prompt de comando, chega-se em uma matriz com duas colunas: a 1 � para k1=1 e a 2 p/ k1=2.
        data.ship(k1).Fmd(:,3) = [Fm6;Fm6(1)];
    end
end


%% Evaluation of retardation functions and infinite added mass

if imemory == 1
    
    A11 =data.hydro.A11;
    A12 =data.hydro.A12;
    A21 =data.hydro.A21;
    A22 =data.hydro.A22;
    B11 =data.hydro.B11;
    B12 =data.hydro.B12;
    B21 =data.hydro.B21;
    B22 =data.hydro.B22;
    
    if ivalida==1
        zeta=.1;
        wn=1;
        for ktime=1:lt
            t = ktime*dt;
            K_teste(1,1,ktime)= wn/((1-zeta^2)^0.5)*exp(-zeta*wn*t)*sin(wn*(1-zeta^2)^0.5*t);
        end
        newdt=1;
        Tij_max = ones(12,1)*1000;
        
        %op��es de excita��o para valida��o
        opcao=1;
        
        switch opcao
            case 1 % funcao degrau unitario
                nu = ones (1,lt);
                %ou
            case 2 %funcao degrau que depois de lt/2 zera
                for cont=1:lt
                    nu(:,cont) = 10;
                    if cont > lt/2
                        nu (:,cont)=0;
                    end
                end
            case 3 %funcao senoide
                %     %ou
                for cont =1:lt
                    nu(:,cont)= sin(pi/5*lt);
                end
        end
        
        for ktime=1:lt
            t = ktime*dt;
            [mu] = convolution_integral(K_teste,Tij_max,ktime,nu,t,newdt);
        end
        
        %curva de compara��o
        % curva para compara��o
        t1=(1:lt)*dt;
        num = [1];
        den = [1 .2 1];
        a=lsim(num,den,nu,t1);
        figure(35)
        i=1;
        plot(variable.t_plot(1:ktime,1),variable.mu(i,1:ktime)','r',variable.t_plot(1:ktime,1),a,'b')
        title('mu (em vermelho) e lsim (em azul), em fun��o do tempo')
        grid on
        %%%%%%
        return
    end
    
    [K11,dK11,Tij11,Nt11,B11,A11] = retardation_function(B11,A11);
    [A11_inf] = infinite_added_mass(A11,K11,dK11);
    
    [K12,dK12,Tij12,Nt12,B12,A12] = retardation_function(B12,A12);
    [A12_inf] = infinite_added_mass(A12,K12,dK12);
    
    [K21,dK21,Tij21,Nt21,B21,A21] = retardation_function(B21,A21);
    [A21_inf] = infinite_added_mass(A21,K21,dK21);
    
    [K22,dK22,Tij22,Nt22,B22,A22] = retardation_function(B22,A22);
    [A22_inf] = infinite_added_mass(A22,K22,dK22);
    
    A_inf = [A11_inf A12_inf ; A21_inf  A22_inf]; % massa adicional total
    
    data.hydro.A_inf = A_inf;
    variable.K11=K11;
    variable.K12=K12;
    variable.K21=K21;
    variable.K22=K22;
    
    data.hydro.Tij11=Tij11;
    data.hydro.Tij12=Tij12;
    data.hydro.Tij21=Tij21;
    data.hydro.Tij22=Tij22;
    
end

%% Time-domain simulation
% initial values
t = 0;

eta0 = [data.ship(1).eta0;data.ship(2).eta0];
nu0 = [data.ship(1).nu0;data.ship(2).nu0];

if icoupling == 1
    sizechitot = data.hydro.sizechi11 + data.hydro.sizechi12 + data.hydro.sizechi21 + data.hydro.sizechi22;
else
    sizechitot = data.hydro.sizechi11 + data.hydro.sizechi22;
end

neq2=18;
nsys = 12 + 12 + sizechitot + 2*0 + neq2; % [eta] [nu] [chi] [control integral terms]
eta_hat0=[data.ship(2).control.xref; data.ship(2).control.yref;data.ship(2).control.psiref];

y = zeros(nsys,lt);
y(:,1) = [eta0;nu0;zeros(neq2,1); zeros(sizechitot,1)]; %aqui o y -> eta ->eta2 ->etaf) %helio
y(25:27,1)=eta_hat0;



variable.ship(1).eta(:,1) = y(1:6,1);
variable.ship(2).eta(:,1) = y(7:12,1);
variable.ship(1).nu(:,1) = y(13:18,1);
variable.ship(2).nu(:,1) = y(19:24,1);
variable.u1p = 0;

% simulation routine

while t <= tfinal
    [t,y(:,ktime+1),nn] = hrkdif(y(:,ktime),dt,t,nn);  %%utiliza��o do hrkdif em que entra o eqsim
    ktime = ktime + 1;
    variable.ship(1).eta(:,ktime) = y(1:6,ktime);
    variable.ship(2).eta(:,ktime) = y(7:12,ktime);
    variable.ship(1).nu(:,ktime) = y(13:18,ktime);
    variable.ship(2).nu(:,ktime) = y(19:24,ktime);
    t
end

texec = cputime - t1;

save tf_graf_d2
% resultados_tf;
figure(1)
plot(tsim,variable.ship(2).waves1st(1,:)/1000)
title('For�a de 1a ordem (PSV) - surge')
xlabel('Tempo (s)')
ylabel('For�a (kN)')
grid on

figure(2)
plot(tsim,variable.ship(2).waves1st(2,:)/1000)
title('For�a de 1a ordem (PSV) - sway')
xlabel('Tempo (s)')
ylabel('For�a (kN)')
grid on

figure(3)
plot(tsim,variable.ship(2).waves1st(6,:)/1000)
title('Momento de 1a ordem (PSV) - yaw')
xlabel('Tempo (s)')
ylabel('Momento (kN.m)')
grid on

figure(4)
plot(tsim,variable.ship(2).tau_wavesmd(1,:)/1000)
title('For�a de 2a ordem (PSV) - surge')
xlabel('Tempo (s)')
ylabel('For�a (kN)')
grid on

figure(5)
plot(tsim,variable.ship(2).tau_wavesmd(2,:)/1000)
title('For�a de 2a ordem (PSV) - sway')
xlabel('Tempo (s)')
ylabel('For�a (kN)')
grid on

figure(6)
plot(tsim,variable.ship(2).tau_wavesmd(6,:)/1000)
title('Momento de 2a ordem (PSV) - yaw')
xlabel('Tempo (s)')
ylabel('Momento (kN.m)')
grid on