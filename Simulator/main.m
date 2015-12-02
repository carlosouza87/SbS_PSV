%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%     Simulation of side-by-side operations between a PSV and an FPSO    %%

clear variables global;close all;clc

global ktime dt data variable
global iwaves1st iwavesmd iwind icontrol_psv ihdsuction isimtype idof icoupling iplot_control icurr imooring
t1 = cputime;
nn=0;

%% Simulation type and number of DOFs
% Determine type of simulation, "isimtype": 
% isimtype = 1, for Cummins equation 
% isimtype = 2, for LF + WF superposition 
isimtype = 1; 

% Determine number o degrees of freedom, "idof":
% idof = 1, for 6 DOF
% idof = 2, for 3 DOF
idof = 1;

%% Simulation data
dimensions      % m-file with definition of ships dimensions
simdata         % Read simulation parameters and organize simulation data into structures

%% Time parameters
if isimtype == 1
    dt = data.constants.dt;
elseif isimtype == 2
    dt = 0.1; % Time step [s]  
    data.constants.dt = dt;
else
    error('Invalid value for isimtype!')
end

tfinal = 100;% total simulation time [s]
lt = tfinal/dt+1;   % number of time steps including zero []
tsim = 0:dt:tfinal; % time vector
ktime = 1;          % time step [s]
variable.dt = dt;
variable.lt = lt;

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

%% Environmental loads
% Waves
load randu1.txt -ascii      % list with uniformly distributed values
spec = 2;                   % flag for wave spectrum (1 = Pierson-Moskowitz, 2 = JONSWAP)
data.ship(1).waves_incid = [data.ship(1).waves_incid;360];
data.ship(2).waves_incid = [data.ship(2).waves_incid;360];

% 1st order wave loads calculation
if iwaves1st == 1
    for k1 = 1:2
        amp = data.ship(k1).w1st_amp;
        pha = data.ship(k1).w1st_pha;
        omg = data.ship(k1).waves_omg;        
        if isimtype == 1
        [Fwf1,Fwf2,Fwf3,Fwf4,Fwf5,Fwf6,Sl1,Sl2,Sl3,Sl4,Sl5,Sl6] = waveloads_1st(tsim,Hs,Tp,spec,amp,pha,omg,randu1);
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
        elseif isimtype == 2
        [Xwf1,Xwf2,Xwf3,Xwf4,Xwf5,Xwf6,Sx1,Sx2,Sx3,Sx4,Sx5,Sx6] = waveloads_1st(tsim,Hs,Tp,spec,amp,pha,omg,randu1);
        data.ship(k1).Xwf(:,:,1) = [Xwf1 Xwf1(:,1)];
        data.ship(k1).Xwf(:,:,2) = [Xwf2 Xwf2(:,1)];
        data.ship(k1).Xwf(:,:,3) = [Xwf3 Xwf3(:,1)];
        data.ship(k1).Xwf(:,:,4) = [Xwf4 Xwf4(:,1)];
        data.ship(k1).Xwf(:,:,5) = [Xwf5 Xwf5(:,1)];
        data.ship(k1).Xwf(:,:,6) = [Xwf6 Xwf6(:,1)];
        data.ship(k1).Sx(:,:,1) = [Sx1 Sx1(:,1)];
        data.ship(k1).Sx(:,:,2) = [Sx2 Sx2(:,1)];
        data.ship(k1).Sx(:,:,3) = [Sx3 Sx3(:,1)];
        data.ship(k1).Sx(:,:,4) = [Sx4 Sx4(:,1)];
        data.ship(k1).Sx(:,:,5) = [Sx5 Sx5(:,1)];
        data.ship(k1).Sx(:,:,6) = [Sx6 Sx6(:,1)];
        end
    end   
end

% 2nd order wave loads calculation
if iwavesmd == 1
    for k1 = 1:2
        amp = data.ship(k1).w2nd_amp;
        pha = data.ship(k1).w2nd_pha;
        omg = data.ship(k1).waves_omg;
        [Fm1,Fm2,Fm6] = meandrift(Hs,Tp,spec,amp,pha,omg);
        data.ship(k1).Fmd(:,1) = [Fm1;Fm1(1)]; 
        data.ship(k1).Fmd(:,2) = [Fm2;Fm2(1)];
        data.ship(k1).Fmd(:,3) = [Fm6;Fm6(1)];
    end
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