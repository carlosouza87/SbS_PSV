% Reads 1st order loads (for isimtype == 1) or motions (for isimtype == 2)
% transfer functions from WAMIT file and save in Nomg X Nincid X 6 
% matrices, where Nomg is the number of frequencies and Nincid the 
% number of incidence considered.

if isimtype == 1
    w1st_file = [caseid '.3']; % Loads are read from WAMIT .3 file
    scl_factor = rho*g; % Scaling factor
    m = 2; % Exponent for ULEN during scaling
elseif isimtype == 2
    w1st_file = [caseid '.4']; % Motions are read from WAMIT .4 file
    scl_factor = 1; % Scaling factor
    m = 0; % Exponent for ULEN during scaling
else
    error('Invalid value for isimtype!')
end

% Read the input file, assigning each column to the corresponding array
[periods,incid,dof,amp,pha,Re,Im] = textread(w1st_file);

periods = unique(periods); % Array of non-repeating periods
incid = unique(incid); % Array of non-repeating incidence directions

[omg,idx] = sort(2*pi./periods);  % "omega" array, for frequency in rad/s
Nomg = length(omg);  % Number of frequencies considered
Nincid = length(incid);   % Number of incidence directions considered
Ndof = length(amp)/(Nomg*Nincid); % Total number of degrees of freedom considered (sum for both ships)

w1st_amp1 = zeros(Nomg,Nincid,6); % Matrix of 1st order force/motion amplitudes for ship 1
w1st_amp2 = zeros(Nomg,Nincid,6); % Matrix of 1st order force/motion amplitudes for ship 2
w1st_pha1 = zeros(Nomg,Nincid,6); % Matrix of 1st order force/motion phases for ship 1
w1st_pha2 = zeros(Nomg,Nincid,6); % Matrix of 1st order force/motion phases for ship 2

% Write the data into Nomg X Nincid X 6 matrices
cont = 0;
for k1 = 1:Nomg
    for k2 = 1:Nincid
        for k3 = 1:6
            % Amplitude for DOFs 1, 2 and 3 are scaled with ULEN^2, for 
            % isimtype== 1, and ULEN^1, for isimtype == 4. For DOFs 4, 5 
            % and 6 the scaling is made with ULEN^3 and ULEN^2, respectively.
            if k3 <= 3                        
                w1st_amp1(Nomg-k1+1,k2,k3) = amp(cont+k3,1)*scl_factor*ULEN^m;
                w1st_amp2(Nomg-k1+1,k2,k3) = amp(cont+k3+6,1)*scl_factor*ULEN^m;            
            else
                w1st_pha1(Nomg-k1+1,k2,k3) = amp(cont+k3,1)*scl_factor*ULEN^(m+1);
                w1st_pha2(Nomg-k1+1,k2,k3) = amp(cont+k3+6,1)*scl_factor*ULEN^(m+1);
            end
            w1st_pha1(Nomg-k1+1,k2,k3) = pha(cont+k3,1);
            w1st_pha2(Nomg-k1+1,k2,k3) = pha(cont+k3+6,1);
        end
        cont = cont + Ndof;
    end
end

% Clear some variables to avoid confusion
clear periods dof amp pha Re Im

% 2nd order loads are read from WAMIT .9 file, which corresponds to the
% direct integration of the pressure over the hull. The data is written in 
% Nomg X Nincid X 6 matrices.
w2nd_file = [caseid '.9']; %  % Loads are read from WAMIT .9 file

% Read the input file, assigning each column to the corresponding array
[periods,incid1,incid2,dof,amp,pha,Re,Im] = textread(w2nd_file);

periods = unique(periods); % Array of non-repeating periods
incid = unique(incid1); % Array of non-repeating incidence directions

[omg,idx] = sort(2*pi./periods);  % "omega" array, for frequency in rad/s
Nomg = length(omg);  % Number of frequencies considered
Nincid = length(incid);   % Number of incidence directions considered
Ndof = length(amp)/(Nomg*Nincid); % Total number of degrees of freedom considered (sum for both ships)

scl_factor = rho*g;  % Scaling factor

w2nd_amp1 = zeros(Nomg,Nincid,6); % Matrix of 2nd order force amplitudes for ship 1
w2nd_amp2 = zeros(Nomg,Nincid,6); % Matrix of 2nd order force amplitudes for ship 2
w2nd_pha1 = zeros(Nomg,Nincid,6); % Matrix of 2nd order phases for ship 1
w2nd_pha2 = zeros(Nomg,Nincid,6); % Matrix of 2nd order phases for ship 2

% Write the data into Nomg X Nincid X 6 matrices
cont = 0;
for k1 = 1:Nomg
    for k2 = 1:Nincid
        for k3 = 1:6
            % Amplitude for DOFs 1, 2 and 3 are scaled with ULEN^1, while 
            % for DOFs 4, 5 and 6 the scaling is made with ULEN^2.
            if k3 <= 3
                w2nd_amp1(Nomg-k1+1,k2,k3) = amp(cont+k3,1)*scl_factor*ULEN;
                w2nd_amp2(Nomg-k1+1,k2,k3) = amp(cont+k3+9,1)*scl_factor*ULEN;
            else
                w2nd_amp1(Nomg-k1+1,k2,k3) = amp(cont+k3,1)*scl_factor*ULEN^2;
                w2nd_amp2(Nomg-k1+1,k2,k3) = amp(cont+k3+9,1)*scl_factor*ULEN^2;
            end
            w2nd_pha1(Nomg-k1+1,k2,k3) = pha(cont+k3,1);
            w2nd_pha2(Nomg-k1+1,k2,k3) = pha(cont+k3+9,1);
        end
        cont = cont + Ndof;
    end
end

% Save relevant variables to be used in the simulations
save('waveloads','w1st_amp1','w1st_amp2','w1st_pha1','w1st_pha2','w2nd_amp1','w2nd_amp2','w2nd_pha1','w2nd_pha2','incid','omg')
