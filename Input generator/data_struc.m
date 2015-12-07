% Saves ships dimensions and hydrodynamic data in a structure named "data"

%% Ships dimensions and properties
% FPSO
% General
data.ship(1).Loa = Loa_1;   % Length overall [m]
data.ship(1).Lpp = Lpp_1;   % Length between perpendiculars [m]
data.ship(1).B = B_1;       % Beam [m]
data.ship(1).T = T_1;       % Draft [m]
data.ship(1).D = D_1;       % Depth [m]
data.ship(1).rg = [0;0;0];   % C.G. position vector [m]
data.ship(1).lL = Cgr_1(1); % Relative distance of CGx to midships [m]
data.ship(1).Cb = Cb_1;      % Block coefficient
data.ship(1).Cy = Cy_1;      % Cross flow coefficient for current sway force calculation
data.ship(1).lCy = lCy_1;    % Cross flow coefficient for current yaw moment calculation
data.ship(1).Vs = Vs_1;   % service speed [m/s]
data.ship(1).S = S_1;         % Wet area [m^2]
data.ship(1).At = At_1;   % Transversal upwater area [m^2]
data.ship(1).Al = Al_1; % Longitudinal upwater area [m^2]
data.ship(1).bow = bow_1;   % Bow type (0 = conventional)
data.ship(1).Mrb = Mrb_1; % Rigid-body inertia matrix
data.ship(1).Ghd = Ghd_1; % Hydrostatic restoration matrix
data.ship(1).GMt = GMt_1; % Transversal metacentric height, FPSO [m]
data.ship(1).GMl = GMl_1; % Longitudinal metacentric height, FPSO [m]
data.ship(1).Cbu = Cbu_1; % Center of buoyancy coordinates [m]
data.ship(1).Cgr = Cgr_1; % Center of gravity coordinates [m]
data.ship(1).Be = Be_1; % External damping matrix

% PSV
% General
data.ship(2).Loa =  Loa_2;   % length overall [m]
data.ship(2).Lpp =  Lpp_2;   % length between perpendiculars [m]
data.ship(2).B = B_2;       % beam [m]
data.ship(2).T = T_2;       % draft [m]
data.ship(2).D = D_2;       % depth [m]
data.ship(2).rg = [0;0;0];   % C.G. position vector [m]
data.ship(2).lL = data.ship(2).rg(2)/data.ship(2).Loa;
data.ship(2).Cb = Cb_2;      % Block coefficient
data.ship(2).Cy = Cy_2;      % Cross flow coefficient for current sway force calculation
data.ship(2).lCy = lCy_2;    % Cross flow coefficient for current yaw moment calculation
data.ship(2).Vs = Vs_2;   % service speed [m/s]
data.ship(2).S = S_2;         % Wet area [m^2]
data.ship(2).At = At_2;   % Transversal upwater area [m^2]
data.ship(2).Al = Al_2; % Longitudinal upwater area [m^2]
data.ship(2).bow = bow_2;   % Bow type (0 = conventional)
data.ship(2).Mrb = Mrb_2; % Rigid-body inertia matrix
data.ship(2).Ghd = Ghd_2; % Hydrostatic restoration matrix
data.ship(2).GMt = GMt_2; % Transversal metacentric height, FPSO [m]
data.ship(2).GMl = GMl_2; % Longitudinal metacentric height, FPSO [m]
data.ship(2).Cbu = Cbu_2; % Center of buoyancy coordinates [m]
data.ship(2).Cgr = Cgr_2; % Center of gravity coordinates [m]
data.ship(2).Be = Be_2; % External damping matrix

% % Controler references and parameters
% data.ship(2).control.Uref = data.ship(2).Vs;% surge reference velocity [m/s]
% data.ship(2).control.yref = 0;% sway reference velocity [m/s]
% data.ship(2).control.Uref_p = 0; % Uref time derivative [m/s] 
% 
% data.ship(2).control.Kp_x = 1e8*1e-1;    % surge speed proportional gain
% data.ship(2).control.Kd_x = 1e3*1e1;    % surge speed derivative gain
% data.ship(2).control.Ki_x = 0*1e7;  % surge speed integral gain
% data.ship(2).control.Kp_y = 1e8;    % sway speed proportional gain
% data.ship(2).control.Kd_y = 1e5;    % sway speed derivative gain
% data.ship(2).control.Ki_y = 0*1e7;  % sway speed integral gain
% data.ship(2).control.Kp_psi = 1e10;       % proportional gain (yaw)
% data.ship(2).control.Kd_psi = 1e10;   % derivative gain (yaw)
% data.ship(2).control.Ki_psi = 0*1e6;      % integral gain (yaw)

%% Hydrodynamic data
load hydro_data

data.constants.dt = dt;
if isimtype == 1
    data.hydro.A11_inf = A11_inf;
    data.hydro.A12_inf = A12_inf;
    data.hydro.A21_inf = A21_inf;
    data.hydro.A22_inf = A22_inf;
    data.hydro.K11 = K11;
    data.hydro.K12 = K12;
    data.hydro.K21 = K21;
    data.hydro.K22 = K22;
elseif isimtype == 2
    data.hydro.A11 = A11;
    data.hydro.A12 = A12;
    data.hydro.A21 = A21;
    data.hydro.A22 = A22;
end

% Wave loads
load waveloads

% Frequency and incidence direction
data.ship(1).waves_omg = omg;
data.ship(1).waves_incid = incid;
data.ship(2).waves_omg = omg;
data.ship(2).waves_incid = incid;

% 1st order
data.ship(1).w1st_amp = w1st_amp1;
data.ship(1).w1st_pha = w1st_pha1;
data.ship(2).w1st_amp = w1st_amp2;
data.ship(2).w1st_pha = w1st_pha2;

% 2nd order
data.ship(1).w2nd_amp = w2nd_amp1;
data.ship(1).w2nd_pha = w2nd_pha1;
data.ship(2).w2nd_amp = w2nd_amp2;
data.ship(2).w2nd_pha = w2nd_pha2;

save (['simdata_' caseid], 'data')