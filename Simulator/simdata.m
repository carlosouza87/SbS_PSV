%% Constant values
knt2ms = 0.5144;  % conversion factor knots - m/s
rho_water = 1025; % water density [kg/m^3]
g = 9.8; % acceleration of gravity [m/s^2]
Vs = 0; % service speed [m/s]

%% WAMIT data 
% Here the user should specify the data to be read, according to the
% loading case, lateral distance between ships and damping factor for the
% lid generalized modes

%LOAD THE DATA FROM INPUT GENERATOR
load hydrodata_cs1d3m_D_3E5
% load hydrodata_fossen

%% Environment data
% wave - 1st order
load randu1.txt -ascii % list with uniformly distributed values
data.environment.betaw = betaw; % waves mean direction [deg]
data.environment.Hs = Hs;   % significant wave height [m]
data.environment.Tp = Tp; % wave modal period [s]
data.environment.w_s = 0:0.05:6;   % frequencies for spectrum calculation

% wind
load dados_OCIMF77.txt -ascii
data.environment.coefwind = dados_OCIMF77; % OCIMF'77 coefficients
data.environment.gammaw = gammaw;    % wind incidence direction [deg]
data.environment.Uw = Uw;    % wind velocity [m/s]
gammaw;
% current
data.environment.alphac = alphac;    % current incidence direction [deg]
data.environment.Uc= Uc;    % current velocity [m/s]


% Fixed frequency matrix 
w = 2*pi/Tp;
for k1 = 1:6
    for k2 = 1:6
        %separando por k1 e k2 e criando um vetor descartável a cada loop
        %para uma faixa de frequência de 1 até esse limite data.hydro.freqs
        for k3 = 1:length(data.hydro.freqs); %k3 representa a variação de frequencia e k1 e k2 a combinação dos graus de liberdade
            A11(k3) = data.hydro.A11(k1,k2,k3);
            A12(k3) = data.hydro.A12(k1,k2,k3);
            A21(k3) = data.hydro.A21(k1,k2,k3);
            A22(k3) = data.hydro.A22(k1,k2,k3);
            B11(k3) = data.hydro.B11(k1,k2,k3);
            B12(k3) = data.hydro.B12(k1,k2,k3);
            B21(k3) = data.hydro.B21(k1,k2,k3);
            B22(k3) = data.hydro.B22(k1,k2,k3);
        end
        %algoritmo para interpolação de uma frequência fixa
        A11_fixfreq(k1,k2) = interp1(data.hydro.freqs,A11,w);%essa frequencia w não é repesentativa!
        A12_fixfreq(k1,k2) = interp1(data.hydro.freqs,A12,w);
        A21_fixfreq(k1,k2) = interp1(data.hydro.freqs,A21,w);
        A22_fixfreq(k1,k2) = interp1(data.hydro.freqs,A22,w);
        B11_fixfreq(k1,k2) = interp1(data.hydro.freqs,B11,w);
        B12_fixfreq(k1,k2) = interp1(data.hydro.freqs,B12,w);
        B21_fixfreq(k1,k2) = interp1(data.hydro.freqs,B21,w);
        B22_fixfreq(k1,k2) = interp1(data.hydro.freqs,B22,w);
    end
end

data.hydro.A11_fixfreq = A11_fixfreq;
data.hydro.A12_fixfreq = A12_fixfreq;
data.hydro.A21_fixfreq = A21_fixfreq;
data.hydro.A22_fixfreq = A22_fixfreq;
data.hydro.B11_fixfreq = B11_fixfreq;
data.hydro.B12_fixfreq = B12_fixfreq;
data.hydro.B21_fixfreq = B21_fixfreq;
data.hydro.B22_fixfreq = B22_fixfreq;


%% Ships data

% FPSO
% General
data.ship(1).Loa = Loa_FPSO;   % length overall [m]
data.ship(1).Lpp = Lpp_FPSO;   % length between perpendiculars [m]
data.ship(1).B = B_FPSO;       % beam [m]
data.ship(1).T = T_FPSO;       % draft [m]
data.ship(1).D = D_FPSO;       % depth [m]
data.ship(1).rg = [0;0;0];   % C.G. position vector [m]
data.ship(1).lL = data.ship(1).rg(2)/data.ship(1).Loa;         %
data.ship(1).Cb = Cb_FPSO;      % block coefficient
data.ship(1).Cy = 0.52;      % cross flow coefficient for current sway force calculation      ?
data.ship(1).lCy = 0.050;     % cross flow coefficient for current yaw moment calculation     ?
% data.ship(1).Vs = 5.0*knt2ms;   % service speed [m/s]
data.ship(1).Vs = Vs;   % service speed [m/s]
data.ship(1).S = S_FPSO;         % wet area [m^2]
data.ship(1).At = At_FPSO;   % transversal upwater area [m^2]
data.ship(1).Al = Al_FPSO; % longitudinal upwater area [m^2]
data.ship(1).bow = 0;   % bow type (0 = conventional)
data.ship(1).load = load_FPSO;  % loading condition (0 = ballast)
data.ship(1).role = 1;  % ship role (1 = guide ship)
dsp = data.ship(1).MRB(1,1)/rho_water;
Awp = Awp_FPSO; % water plane area (full loaded) ship [m]
GMt = GMt_FPSO; % transversal metacentric height (ballasted ship) [m]%considering KG in the water line
GMl = GMl_FPSO; % longitudinal metacentric height (ballasted ship) [m]
data.ship(1).G =diag([0;0;rho_water*g*Awp;rho_water*g*dsp*GMt;rho_water*g*dsp*GMl;0]);

% Initial state vectors
data.ship(1).eta0 = [0;0;0;0*pi/180;0;0];
data.ship(1).nu0 = [data.ship(1).Vs;0;0;0;0;0];

% Control
data.ship(1).control.Uref = data.ship(1).Vs;% surge reference velocity [m/s]
data.ship(1).control.yref = 0;% sway reference velocity [m/s]
data.ship(1).control.Uref_p = 0; % Uref time derivative [m/s]


data.ship(1).control.Kp_x = 1e8*1e-1;    % surge speed proportional gain
data.ship(1).control.Kd_x = 1e3*1e1;    % surge speed derivative gain
data.ship(1).control.Ki_x = 0*1e7;  % surge speed integral gain
data.ship(1).control.Kp_y = 1e8;    % sway speed proportional gain
data.ship(1).control.Kd_y = 1e5;    % sway speed derivative gain
data.ship(1).control.Ki_y = 0*1e7;  % sway speed integral gain
data.ship(1).control.Kp_psi = 1e10;       % proportional gain (yaw)
data.ship(1).control.Kd_psi = 1e10;   % derivative gain (yaw)
data.ship(1).control.Ki_psi = 0*1e6;      % integral gain (yaw)

% Parameters for mean drift forces calculation scheduling
variable.ship(1).calc_mdrift = 1;
variable.ship(1).u_md = data.ship(1).nu0(1) + 1;
variable.ship(1).v_md = data.ship(1).nu0(2) + 1;
variable.ship(1).psi_md = data.ship(1).eta0(6) + 10*pi/180;

% Parameters for wind forces calculation scheduling
variable.ship(1).calc_wind = 1;
variable.ship(1).u_wn = data.ship(1).nu0(1) + 1;
variable.ship(1).v_wn = data.ship(1).nu0(2) + 1;
variable.ship(1).psi_wn = data.ship(1).eta0(6) + 10*pi/180;

% PSV
% General
data.ship(2).Loa =  Loa_PSV;   % length overall [m]
data.ship(2).Lpp =  Lpp_PSV;   % length between perpendiculars [m]
data.ship(2).B = B_PSV;       % beam [m]
data.ship(2).T = T_PSV;       % draft [m]
data.ship(2).D = D_PSV;       % depth [m]
data.ship(2).rg = [0;0;0];   % C.G. position vector [m]
data.ship(2).lL = data.ship(2).rg(2)/data.ship(2).Loa;
data.ship(2).Cb = Cb_PSV;      % block coefficient
data.ship(2).Cy = 0.70;      % cross flow coefficient for current sway force calculation       
data.ship(2).lCy = 0.048;     % cross flow coefficient for current yaw moment calculation      
% data.ship(2).Vs = 5.0*knt2ms;   % service speed [m/s]
data.ship(2).Vs = Vs; %2.5720 m/s   % service speed [m/s]
data.ship(2).S = S_PSV;            % wet area [m^2]
data.ship(2).At = At_PSV; % Transversal upwater area [m^2]
data.ship(2).Al = Al_PSV; % Longitudinal upwater area [m^2]
data.ship(2).bow = 1; % Bow type (1 = bulbous)                                                
data.ship(2).load = load_PSV; % Loading condition (1 = full)
data.ship(2).role = 0; % ship role (0 = towed ship)
dsp = data.ship(2).MRB(1,1)/rho_water;
Awp = Awp_PSV; % Water plane area (full loaded) ship [m]
GMt = GMt_PSV; % Transversal metacentric height (full loaded ship) [m]
GMl = GMl_PSV; % Longitudinal metacentric height (full loaded ship) [m]
data.ship(2).G = diag([0;0;rho_water*g*Awp;rho_water*g*dsp*GMt;rho_water*g*dsp*GMl;0]);


% Initial state vectors
 data.ship(2).eta0 = [0;46.75;0;0;0;0]; %CG foi colocado na linha d'água, foi trocado de 55 para 40, soma das duas meia bocas mais distancia entre navios
 data.ship(2).nu0 = [data.ship(2).Vs;0;0;0;0;0];

% Parameters for mean drift forces calculation scheduling
variable.ship(2).calc_mdrift = 1;
variable.ship(2).u_md = data.ship(2).nu0(1) + 1;
variable.ship(2).v_md = data.ship(2).nu0(2) + 1;
variable.ship(2).psi_md = data.ship(2).eta0(6) + 10*pi/180;

% Parameters for wind forces calculation scheduling
variable.ship(2).calc_wind = 1;
variable.ship(2).u_wn = data.ship(2).nu0(1) + 1;
variable.ship(2).v_wn = data.ship(2).nu0(2) + 1;
variable.ship(2).psi_wn = data.ship(2).eta0(6) + 10*pi/180;


data.ship(2).control.xref=0;
data.ship(2).control.yref=46.75;
data.ship(2).control.psiref=0;

data.ship(2).control.Kp_x = 1e8*1e-2*.5/10;%*1e-1;%*1e-2;    % surge speed proportional gain
data.ship(2).control.Kd_x = 1e6/10;%1e3*1e-1*1e2*1e-3;    % surge speed derivative gain
data.ship(2).control.Ki_x = 0*1e7;  % surge speed integral gain

data.ship(2).control.Kp_y = 1e8*1e-2*.5/10;%*2;    % sway speed proportional gain
data.ship(2).control.Kd_y = 1e5*1e-2*2*10*10/10;    % sway speed derivative gain
data.ship(2).control.Ki_y = 0*1e7;  % sway speed integral gain

data.ship(2).control.Kp_psi = 1e10*1e-2;       % proportional gain (yaw)
data.ship(2).control.Kd_psi = 1e10*1e-2;   % derivative gain (yaw)
data.ship(2).control.Ki_psi = 0*1e6;      % integral gain (yaw)


