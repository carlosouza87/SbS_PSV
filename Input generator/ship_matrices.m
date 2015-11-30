% CONTINUAR DE ROLL DAMPING!

% Define rigid-body inertia matrices for both vessels
MRB1 = [318620 0 0 0 0 0;
        0 318620 0 0 0 0;
        0 0 318620 0 0 0;
        0 0 0 1.20e8 0 0;
        0 0 0 0 1.91e9 0;
        0 0 0 0 0 1.91e9]*1e3; % FPSO

MRB2 = [8170 0 0 0 0 0;
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

data.ship(1).MRB = MRB1;
data.ship(2).MRB = MRB2;


% Define hydrostatic restoration matrices for both vessels
scl = rho*g;
L = ULEN;

C1 = [0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 15705*scl*L^2 0 0.14358e6*scl*L^3 0;
      0 0 0 0.24356e7*scl*L^4 0 -0.10537e6*scl*L^4;
      0 0 0.14358e6*scl*L^3 0 0.11347e9*scl*L^4 0;
      0 0 0 0 0 0]; % FPSO
 
C2 = [0 0 0 0 0 0; 
      0 0 0 0 0 0;
      0 0 1458.7*scl*L^2 0 7177.1*scl*L^3 0;
      0 0 0 24197*scl*L^4 0 228.09*scl*L^4;
      0 0 7177.1*scl*L^3 0 0.75337e6*scl*L^4 0;
      0 0 0 0 0 0]; % PSV
  
  data.ship(1).G = C1;
  data.ship(2).G = C2;
      
      
% External roll damping


%--------------------------------------------------------------------------
% Read rigid-body mass parameters from *.out file
%--------------------------------------------------------------------------
count = 0;

%fid1 = fopen(strcat(vessel1,'.out'));
fid1 = fopen(strcat('FPSO.out'));

rho = 1025;
g = 9.8;


while feof(fid1) == 0,
    count = count + 1;
    txt = char(fgetl(fid1));
    
    % GDF file length scale
    if  strfind(txt,'Length scale:')
        ULEN_1 = str2num(txt(52:length(txt)));
    end
    
    if  strfind(txt,' C(3,3),C(3,4),C(3,5):')
        
        C3_1 = str2num(txt(23:length(txt)));
        txt = char(fgetl(fid1));
        C4_1 = str2num(txt(23:length(txt)));
        txt = char(fgetl(fid1));
        C5_1 = str2num(txt(23:length(txt)));
        
        % scaling to SI units (Wamit manual p. 4-2)        
        C3_1 = C3_1 .* rho*g .* [ULEN_1^2 ULEN_1^3 ULEN_1^3];
        C4_1 = C4_1 .* rho*g .* [ULEN_1^4 ULEN_1^4 ULEN_1^4];
        C5_1 = C5_1 .* rho*g .* [ULEN_1^4 ULEN_1^4];
        
        % spring stiffness matrix in global coordinates (Wamit axes)
        % Wamit manual p. 4-2
        C_wamit_1 = zeros(6,6);
        C_wamit_1(3:6,3:6) =...
            [ C3_1 0
            C3_1(2) C4_1
            C3_1(3) C4_1(2) C5_1
            0 0 0 0 ];
        
        % spring stiffness matrix in CO (Fossen axes)
        C_wamit_1 = Tscale*C_wamit_1*Tscale;
        data.ship(1).G = C_wamit_1;
        data.ship(1).GM_T  = C_wamit_1(4,4)/(MRB1(1)*g);
        data.ship(1).GM_L  = C_wamit_1(5,5)/(MRB1(1)*g);
        
    end
    
    % CG and CB
    if  strfind(txt,'Center of Buoyancy')
        temp = str2num(txt(33:length(txt)));
        C_B_1 = T_gdf.*temp;   % Fossen axes
        data.ship(1).CB = [C_B_1(1) C_B_1(2) draft1-C_B_1(3)];
    end
    if  strfind(txt,'Center of Gravity')
        temp = str2num(txt(33:length(txt)));
        C_G_1 = T_gdf.*temp;   % Fossen axes
        data.ship(1).CG = [C_G_1(1) C_G_1(2) draft1-C_G_1(3)];
    end
    
end   % End WHILE

fclose(fid1);

count = 0;

%fid2 = fopen(strcat(vessel2,'.out'));
fid2 = fopen(strcat('PSV.out'));

while feof(fid2) == 0,
    count = count + 1;
    txt = char(fgetl(fid2));
    
    % GDF file length scale
    if  strfind(txt,'Length scale:')
        ULEN_2 = str2num(txt(52:length(txt)));
    end
    
    if  strfind(txt,' C(3,3),C(3,4),C(3,5):')
        
        C3_2 = str2num(txt(23:length(txt)));
        txt = char(fgetl(fid2));
        C4_2 = str2num(txt(23:length(txt)));
        txt = char(fgetl(fid2));
        C5_2 = str2num(txt(23:length(txt)));
        
        % scaling to SI units (Wamit manual p. 4-2)        
        C3_2 = C3_2 .* rho*g .* [ULEN_2^2 ULEN_2^3 ULEN_2^3];
        C4_2 = C4_2 .* rho*g .* [ULEN_2^4 ULEN_2^4 ULEN_2^4];
        C5_2 = C5_2 .* rho*g .* [ULEN_2^4 ULEN_2^4];
        
        % spring stiffness matrix in global coordinates (Wamit axes)
        % Wamit manual p. 4-2
        C_wamit_2 = zeros(6,6);
        C_wamit_2(3:6,3:6) =...
            [ C3_2 0
            C3_2(2) C4_2
            C3_2(3) C4_2(2) C5_2
            0 0 0 0 ];
        
        % spring stiffness matrix in CO (Fossen axes)
        C_wamit_2 = Tscale*C_wamit_2*Tscale;        
        data.ship(2).G = C_wamit_2;
        data.ship(2).GM_T  = C_wamit_2(4,4)/(MRB1(1)*g);
        data.ship(2).GM_L  = C_wamit_2(5,5)/(MRB1(1)*g);        
    end
    
    % CG and CB
    if  strfind(txt,'Center of Buoyancy')
        temp = str2num(txt(33:length(txt)));
        C_B_2 = T_gdf.*temp;   % Fossen axes
        data.ship(2).CB = [C_B_2(1) C_B_2(2) draft1-C_B_2(3)];
    end
    if  strfind(txt,'Center of Gravity')
        temp = str2num(txt(33:length(txt)));
        C_G_2 = T_gdf.*temp;   % Fossen axes
        data.ship(2).CG = [C_G_2(1) C_G_2(2) draft1-C_G_2(3)];
    end
    
end   % End WHILE

fclose(fid2);

%% Hydrodynamic data
load hydro_data
% load tankerABC

cpl_mem = [1,1;1,3;1,5;2,2;2,4;2,6;3,1;3,3;3,5;4,2;4,4;4,6;5,1;5,3;5,5;6,2;6,4;6,6];

data.hydro.freqs = freqs;
data.hydro.A11 = A11;
data.hydro.A12 = A12;
data.hydro.A21 = A21;
data.hydro.A22 = A22;
data.hydro.B11 = B11;
data.hydro.B12 = B12;
data.hydro.B21 = B21;
data.hydro.B22 = B22;
data.hydro.sizechi11 = size(data.hydro.Amem11,1);
data.hydro.sizechi12 = size(data.hydro.Amem12,1);
data.hydro.sizechi21 = size(data.hydro.Amem21,1);
data.hydro.sizechi22 = size(data.hydro.Amem22,1);
data.hydro.dof_unst11 = dof_unst11;
data.hydro.dof_unst12 = dof_unst12;
data.hydro.dof_unst21 = dof_unst21;
data.hydro.dof_unst22 = dof_unst22;

%% Wave forces
load waveforces

% 1st order
data.ship(1).FTF_amp = FTF_amp1;
data.ship(1).FTF_pha = FTF_pha1;
data.ship(2).FTF_amp = FTF_amp2;
data.ship(2).FTF_pha = FTF_pha2;
data.ship(1).wavefreqs = freqs;
data.ship(1).waveincid = incid;
data.ship(2).wavefreqs = freqs;
data.ship(2).waveincid = incid;

% 2nd order
data.ship(1).drift_amp = drift_amp1;
data.ship(1).drift_pha = drift_pha1;
data.ship(2).drift_amp = drift_amp2;
data.ship(2).drift_pha = drift_pha2;

save (filename, 'data')