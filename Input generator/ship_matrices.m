% The user determines the matrices of rigid-body inertia, hydrostatic 
% coefficients, and external linear damping, and the coordinates of the centers 
% of buoyancy and gravity. Then, the program loads the previously processed data
% for radiation and wave loads, and saves everything in a structure named "ship",
% which will be used in the simulations.

% Define rigid-body inertia matrices for both vessels
Mrb1 = [318620 0 0 0 0 0;
        0 318620 0 0 0 0;
        0 0 318620 0 0 0;
        0 0 0 1.20e8 0 0;
        0 0 0 0 1.91e9 0;
        0 0 0 0 0 1.91e9]*1e3; % FPSO

Mrb2 = [8170 0 0 0 0 0;
        0 8170 0 0 0 0;
        0 0 8170 0 0 0;
        0 0 0 3.62e5 0 0;
        0 0 0 0 4.04e6 0;
        0 0 0 0 0 4.04e6]*1e3; % PSV

T_gdf = [1 -1 -1];             % Wamit2Fossen axes (sign correction)
Tscale = diag([T_gdf T_gdf]);  % 6 DOF transformation matrix for A and B data


% % Transform inertia matrices to Fossen axes
% MRB1 =  Tscale*MRB1*Tscale;
% MRB2 =  Tscale*MRB2*Tscale;

data.ship(1).Mrb= Mrb1;
data.ship(2).Mrb = Mrb2;


% Define hydrostatic restoration matrices for both vessels
scl = rho*g;
L = ULEN;

Ghd1 = [0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 15705*scl*L^2 0 0.14358e6*scl*L^3 0;
      0 0 0 0.24356e7*scl*L^4 0 -0.10537e6*scl*L^4;
      0 0 0.14358e6*scl*L^3 0 0.11347e9*scl*L^4 0;
      0 0 0 0 0 0]; % FPSO
 
Ghd2 = [0 0 0 0 0 0; 
      0 0 0 0 0 0;
      0 0 1458.7*scl*L^2 0 7177.1*scl*L^3 0;
      0 0 0 24197*scl*L^4 0 228.09*scl*L^4;
      0 0 7177.1*scl*L^3 0 0.75337e6*scl*L^4 0;
      0 0 0 0 0 0]; % PSV
  
data.ship(1).Ghd = Ghd1;  
data.ship(2).Ghd = Ghd2;

% Centers of buoyancy and gravity
Cbu1 = [0.34;0;-3.34]; % Center of buoyancy, FPSO [m]
Cbu2 = [-0.03;0;-2.00]; % Center of buoyancy, PSV [m]

Cgr1 = [0;0;0];  % Center of gravity, FPSO [m]
Cgr2 = [0;0;0]; % Center of gravity, PSV [m]

data.ship(1).Cbu = Cbu1;
data.ship(1).Cgr = Cgr1;
data.ship(2).Cbu = Cbu2;
data.ship(2).Cgr = Cgr2;

% Metacentric heights
data.ship(1).GM_T = Ghd1(4,4)/(Mrb1(1)*g); % Transversal metacentric height, FPSO [m]
data.ship(1).GM_L = Ghd1(5,5)/(Mrb1(1)*g); % Longitudinal metacentric height, FPSO [m]

data.ship(2).GM_T = Ghd2(4,4)/(Mrb2(1)*g); % Transversal metacentric height, PSV [m]
data.ship(2).GM_L = Ghd2(5,5)/(Mrb2(1)*g); % Longitudinal metacentric height, PSV [m]
      
% External roll damping
Be1 = zeros(6,6);
Be1(4,4) = 6.0971E+06; % Only external damping in roll is considered [N.s]
data.ship(1).Be = Be1;

Be2 = zeros(6,6);
Be2(4,4) = 2.4635E+04; % Only external damping in roll is considered [N.s]
data.ship(2).Be = Be2;


%% Hydrodynamic data
load hydro_data

if imemory == 0
data.hydro.A11 = A11;
data.hydro.A12 = A12;
data.hydro.A21 = A21;
data.hydro.A22 = A22;
elseif imemory == 1
data.hydro.A11_inf = A11_inf;
data.hydro.A12_inf = A12_inf;
data.hydro.A21_inf = A21_inf;
data.hydro.A22_inf = A22_inf;
data.hydro.K11 = K11;
data.hydro.K12 = K12;
data.hydro.K21 = K21;
data.hydro.K22 = K22;
end

%% Wave forces
load waveforces

% 1st order
data.ship(1).w1st_amp = w1st_amp1;
data.ship(1).w1st_pha = w1st_pha1;
data.ship(2).w1st_amp = w1st_amp2;
data.ship(2).w1st_pha = w1st_pha2;
% CHECAR WAVELOADS!
data.ship(1).wavefreqs = freqs;
data.ship(1).waveincid = incid;
data.ship(2).wavefreqs = freqs;
data.ship(2).waveincid = incid;

% 2nd order
data.ship(1).w2nd_amp = w2nd_amp1;
data.ship(1).w2nd_pha = w2nd_pha1;
data.ship(2).w2nd_amp = w2nd_amp2;
data.ship(2).w2nd_pha = w2nd_pha2;

save (caseid, 'data')