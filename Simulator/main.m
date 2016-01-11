%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%     Simulation of side-by-side operations between a PSV and an FPSO    %%

clear variables global;close all;clc

global ktime dt data variable
global flag
t1 = cputime;
nn=0;

%% Input caseid
caseid = 'conjunto_10';

%% Time parameters
dt = 0.1;% time step [s]   
tfinal = 600;% total simulation time [s]
lt = tfinal/dt;   % number of time steps including zero []
tsim = 0:dt:tfinal; % time vector
ktime = 1;          % time step [s]
data.constants.dt = dt;
variable.dt = dt;
variable.lt=lt;

%% Simulation type and number of DOFs
% Determine type of simulation, "isimtype": 
% isimtype = 1, for Cummins equation 
% isimtype = 2, for LF + WF superposition 
flag.isimtype = 2;   

% Determine number o degrees of freedom, "idof":
% idof = 1, for 6 DOF
% idof = 2, for 3 DOF
flag.idof = 1;

%% Flags for switching modules on/off
flag.iwaves1st = 1;       % Flag for 1st order wave loads
flag.iwavesmd = 0;        % Flag for meand drift loads
flag.iwind = 0;           % Flag for wind loads
flag.icurr = 0;           % Flag for current loads
flag.icontrol_psv = 0;    % Flag for the psv control system
flag.ihdsuction = 0;      % Flag for hydrodynamci suction loads
flag.imooring = 0;        % Flag for FPSO mooring system
% flag.icoupling = 0;       % Flag for coupled dynamics
% imemory = 0;       % Flag for memory effects
% ivalida = 0;       % Flag to confirm validation of convolution evaluation


%% Loading of simulation data
simdata         % Read simulation parameters and organize simulation data into structures

%% Environmental loads
% Waves

data.ship(1).waves_incid = [data.ship(1).waves_incid;360];
data.ship(2).waves_incid = [data.ship(2).waves_incid;360];

% 1st order wave loads calculation
if flag.iwaves1st == 1
    for k1 = 1:2
        amp = data.ship(k1).w1st_amp;
        pha = data.ship(k1).w1st_pha;
        freqs = data.ship(k1).waves_omg;
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
if flag.iwavesmd == 1
    for k1 = 1:2
        amp = data.ship(k1).w2nd_amp;
        pha = data.ship(k1).w2nd_pha;
        freqs = data.ship(k1).waves_omg;
        [Fm1,Fm2,Fm6] = meandrift(Hs,Tp,spec,amp,pha,freqs);
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

% Ensure that the integration vector size matches the size of all vectors
% to be integrated in eqsim.
nsys = 12 + 12; % [eta] [nu]

y = zeros(nsys,lt);
y(:,1) = [eta0;nu0];

variable.ship(1).eta(:,1) = y(1:6,1);
variable.ship(2).eta(:,1) = y(7:12,1);
variable.ship(1).nu(:,1) = y(13:18,1);
variable.ship(2).nu(:,1) = y(19:24,1);

% simulation routine

while t <= tfinal
    [t,y(:,ktime+1)] = hrkdif(y(:,ktime),dt,t);  
    ktime = ktime + 1;
    variable.ship(1).eta(:,ktime) = y(1:6,ktime);
    variable.ship(2).eta(:,ktime) = y(7:12,ktime);
    variable.ship(1).nu(:,ktime) = y(13:18,ktime);
    variable.ship(2).nu(:,ktime) = y(19:24,ktime);
    t
end

if flag.isimtype == 1
    save mem_on variable
else
    save sup variable
end

texec = cputime - t1;

