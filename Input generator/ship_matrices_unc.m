% clear variables;close all;clc
%
% vessel1  = input('Vessel 1 file name: ','s');
% draft1  = input('Vessel 1 draft: ');
% vessel2  = input('Vessel 2 file name: ','s');
% draft2  = input('Vessel 2 draft: ');

eval(['load MRB_' vessel1 '.txt -ascii']);
eval(['MRB1 = MRB_' vessel1 ' ;'])
eval(['load MRB_' vessel2 '.txt -ascii'])
eval(['MRB2 = MRB_' vessel2 ';'])

T_gdf = [1 -1 -1];             % Wamit2Fossen axes (sign correction)
Tscale = diag([T_gdf T_gdf]);  % 6 DOF transformation matrix for A and B data


%% Rigid body data
% mass matrices in CO (Fossen axes)
MRB1 =  Tscale*MRB1*Tscale;
MRB2 =  Tscale*MRB2*Tscale;
data.ship(1).MRB = MRB1;
data.ship(2).MRB = MRB2;

%--------------------------------------------------------------------------
%% Read rigid-body mass parameters from *.out file
%--------------------------------------------------------------------------
count = 0;

fid1 = fopen(strcat(vessel1,'.out'));

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
    
    % Infinite-frequency added mass
    if strfind(txt,'Wave period = zero')
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        for k1_smu = 1:18
            a = str2num(char(fgetl(fid1)));
            A11inf(a(1),a(2)) = rho*a(3);
        end
        a = 1;
    end
    
    % Zero-frequency added mass
    if strfind(txt,'Wave period = infinite')
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        for k1_smu = 1:18
            a = str2num(char(fgetl(fid1)));
            A11zero(a(1),a(2)) = rho*a(3);
        end
        a = 1;
    end
    
end   % End WHILE

fclose(fid1);

count = 0;

fid2 = fopen(strcat(vessel2,'.out'));

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
    
        % Infinite-frequency added mass
    if strfind(txt,'Wave period = zero')
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        for k1_smu = 1:18
            a = str2num(char(fgetl(fid1)));
            A22inf(a(1),a(2)) = rho*a(3);
        end
        a = 1;
    end
    
    % Zero-frequency added mass
    if strfind(txt,'Wave period = infinite')
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        fgetl(fid1);
        for k1_smu = 1:18
            a = str2num(char(fgetl(fid1)));
            A22zero(a(1),a(2)) = rho*a(3);
        end
        a = 1;
    end
    
end   % End WHILE

fclose(fid2);

%% Hydrodynamic data
load memory_ss
% load tankerABC
%
% A11r = vesselABC.Ar;
% B11r = vesselABC.Br;
% C11r = vesselABC.Cr;
cpl_mem = [1,1;2,2;2,4;2,6;3,3;3,5;4,2;4,4;4,6;5,1;5,5;6,2;6,4;6,6];
[Amem11,Bmem11,Cmem11,sizechi11,dof_mem11] = mem_matrices(cpl_mem,A11r,B11r,C11r);
% [Amem12,Bmem12,Cmem12,sizechi12,dof_mem12] = mem_matrices(cpl_mem,A12r,B12r,C12r);
% [Amem21,Bmem21,Cmem21,sizechi21,dof_mem21] = mem_matrices(cpl_mem,A21r,B21r,C21r);

% A22r = vesselABC.Ar;
% B22r = vesselABC.Br;
% C22r = vesselABC.Cr;
[Amem22,Bmem22,Cmem22,sizechi22,dof_mem22] = mem_matrices(cpl_mem,A22r,B22r,C22r);

data.hydro.Amem11 = Amem11;
data.hydro.Bmem11 = Bmem11;
data.hydro.Cmem11 = Cmem11;
data.hydro.dof_mem11 = dof_mem11;
% data.hydro.Amem12 = Amem12;
% data.hydro.Bmem12 = Bmem12;
% data.hydro.Cmem12 = Cmem12;
% data.hydro.dof_mem12 = dof_mem12;
% data.hydro.Amem21 = Amem21;
% data.hydro.Bmem21 = Bmem21;
% data.hydro.Cmem21 = Cmem21;
% data.hydro.dof_mem21 = dof_mem21;
data.hydro.Amem22 = Amem22;
data.hydro.Bmem22 = Bmem22;
data.hydro.Cmem22 = Cmem22;
data.hydro.dof_mem22 = dof_mem22;
data.hydro.freqs = freqs;
data.hydro.A11(:,:,1) = A11zero;
for k1 = 1:16
    data.hydro.A11(:,:,k1+1) = A11(:,:,k1);
end
data.hydro.A11(:,:,18) = A11inf;
% data.hydro.A12 = A12;
% data.hydro.A21 = A21;
data.hydro.A22(:,:,1) = A22zero;
for k1 = 1:12
    data.hydro.A22(:,:,k1+1) = A22(:,:,k1);
end
data.hydro.A22(:,:,18) = A22inf;
data.hydro.B11(:,:,1) = zeros(6);
for k1 = 1:16
    data.hydro.B11(:,:,k1+1) = B11(:,:,k1);
end
data.hydro.B11(:,:,18) = zeros(6);
% data.hydro.B12 = B12;
% data.hydro.B21 = B21;
data.hydro.B22(:,:,1) = zeros(6);
for k1 = 1:16
    data.hydro.B22(:,:,k1+1) = B22(:,:,k1);
end
data.hydro.B22(:,:,18) = zeros(6);
data.hydro.sizechi11 = size(data.hydro.Amem11,1);
% data.hydro.sizechi12 = size(data.hydro.Amem12,1);
% data.hydro.sizechi21 = size(data.hydro.Amem21,1);
data.hydro.sizechi22 = size(data.hydro.Amem22,1);
data.hydro.dof_unst11 = dof_unst11;
% data.hydro.dof_unst12 = dof_unst12;
% data.hydro.dof_unst21 = dof_unst21;
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

% save coupled data
%
% clear all