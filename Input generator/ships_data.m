% Input dimensions, matrices of rigid-body inertia, hydrostatic
% coefficients, and external linear damping, and the coordinates of the 
% centers of buoyancy and gravity. Then, the program loads the previously 
% processed data for radiation and wave loads, and saves everything in a 
% structure named "ship", which will be used in the simulations.

%% FPSO (full)
Loa_1 = 333.0;   % Length overall [m]
Lpp_1 = 320.0;   % Length between perpendiculars [m]
B_1 = 54.5;      % Beam [m]
T_1 = 21.5;      % Draft [m]
D_1 = 27.8;      % Depth [m]
Delta_1 = 318620;    % Displacement [t]
Nabla_1 = 310850;    % Volumetric displacement [m^3]
Cb_1 = 0.85;     % Block coefficient []
Cy_1 = 0.52;     % Cross flow coefficient for current sway force calculation
lCy_1 = 0.050;   % Cross flow coefficient for current yaw moment calculation
S_1 = 27447;     % Wet area [m^2]
Awp_1 = 15710;   % Waterline area [m^2]
Vs_1 = 0;        % Service speed [m/s]
% GMt_1 = 7.34;    % Transversal metacentric height [m]
% GMl_1 = 370.7;   % Longitudinal metacentric height [m]
At_1 = 3939.0;   % Wind frontal area [m^2]
Al_1 = 14368.0;  % Wind lateral area [m^2]
bow_1 = 0;       % Bow type (0 = conventional)
load_1 = 1;      % Loading condition [1 = full]
Cbu_1 = [0.34;0;-3.34]; % Center of buoyancy [m]
Cgr_1 = [0;0;0];  % Center of gravity [m]

% Rigid-body inertia matrix (usually available in the WAMIT .out file)
Mrb_1 = [318620 0 0 0 0 0;
    0 318620 0 0 0 0;
    0 0 318620 0 0 0;
    0 0 0 1.20e8 0 0;
    0 0 0 0 1.91e9 0;
    0 0 0 0 0 1.91e9]*1e3;

% % Transform inertia matrix to Fossen axes
% T_gdf = [1 -1 -1];             % Wamit2Fossen axes (sign correction)
% Tscale = diag([T_gdf T_gdf]);  % 6 DOF transformation matrix for A and B data
% Mrb1 =  Tscale*MRB1*Tscale;

% Hydrostatic restoration matrix (usually available in the WAMIT .out file)
scl = rho*g;
L = ULEN;

Ghd_1 = [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 15705*scl*L^2 0 0.14358e6*scl*L^3 0;
    0 0 0 0.24356e7*scl*L^4 0 -0.10537e6*scl*L^4;
    0 0 0.14358e6*scl*L^3 0 0.11347e9*scl*L^4 0;
    0 0 0 0 0 0]; 
% Metacentric heights based on hidrostatic coefficients
GMt_1 = Ghd_1(4,4)/(Mrb_1(1)*g); % Transversal metacentric height [m]
GMl_1 = Ghd_1(5,5)/(Mrb_1(1)*g); % Longitudinal metacentric height [m]

% External roll damping
Be_1 = zeros(6,6);
Be_1(4,4) = 6.0971E+06; % Only external damping in roll is considered [N.s]

%% PSV (full)
Loa_2 = 88.8;     % Length overall [m]
Lpp_2 = 80.0;     % Length between perpendiculars [m]
B_2 = 19.0;       % Beam [m]
T_2 = 6.6;        % Draft [m]
D_2 = 8.0;        % Depth [m]
Delta_2 = 8193.8; % Displacement [t]
Nabla_2 = 7994;   % Volumetric displacement [m^3]
Cb_2 = 0.73;      % Block coefficient []
Cy_2 = 0.52;     % Cross flow coefficient for current sway force calculation
lCy_2 = 0.050;   % Cross flow coefficient for current yaw moment calculation
S_2 = 2341.7;     % Wet area [m^2]
Awp_2 = 1454.5;   % Waterline area [m^2]
Vs_2 = 0;        % Service speed [m/s]
% GMt_2 = 3.14;     % Transversal metacentric height [m]
% GMl_2 = 90.9;     % Longitudinal metacentric height [m]
At_2 = 238.45;    % Wind frontal area [m^2]
Al_2 = 339.2;     % Wind lateral area [m^2]
bow_2 = 0;       % Bow type (0 = conventional)
load_2 = 1;       % Loading condition [1 = full]
Cbu_2 = [-0.03;0;-2.00]; % Center of buoyancy, PSV [m]
Cgr_2 = [0;0;0]; % Center of gravity, PSV [m]

% Rigid-body inertia matrix (usually available at the WAMIT .out file)
Mrb_2 = [8170 0 0 0 0 0;
    0 8170 0 0 0 0;
    0 0 8170 0 0 0;
    0 0 0 3.62e5 0 0;
    0 0 0 0 4.04e6 0;
    0 0 0 0 0 4.04e6]*1e3;

% % Transform inertia matrices to Fossen axes
% T_gdf = [1 -1 -1];             % Wamit2Fossen axes (sign correction)
% Tscale = diag([T_gdf T_gdf]);  % 6 DOF transformation matrix for A and B data
% Mrb1 =  Tscale*MRB1*Tscale;

% Hydrostatic restoration matrix (usually available in the WAMIT .out file)
scl = rho*g;
L = ULEN;

Ghd_2 = [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 1458.7*scl*L^2 0 7177.1*scl*L^3 0;
    0 0 0 24197*scl*L^4 0 228.09*scl*L^4;
    0 0 7177.1*scl*L^3 0 0.75337e6*scl*L^4 0;
    0 0 0 0 0 0]; 

% Metacentric heights based on hidrostatic coefficients
GMt_2 = Ghd_2(4,4)/(Mrb_2(1)*g); % Transversal metacentric height, PSV [m]
GMl_2 = Ghd_2(5,5)/(Mrb_2(1)*g); % Longitudinal metacentric height, PSV [m]

% External roll damping
Be_2 = zeros(6,6);
Be_2(4,4) = 2.4635E+04; % Only external damping in roll is considered [N.s]
