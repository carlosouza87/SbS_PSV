clear all;close all;clc

% It is admitted that the WAMIT files bring either both asymptotic values (omg
% = 0 rad/s and omg = Inf) or none of them.

caseid = 'conjunto';
inpt1 = [caseid '.1']; % Input file
imemory = 1;

% Read the input file, assigning each column to the corresponding array
[periods,dof_i,dof_j,A,B] = textread(inpt1);

Nperiods = length(periods); % Number of periods, considering repeated
unique_periods = unique(periods); % Array of non-repeating periods

% Prepare matrices for scaling of added mass and damping matrices (WAMIT
% v 6.2 manual, p. 4-3)
% Aij = Aij' * rho * ULEN^k
% Bij = Bij' * rho * w * ULEN^k
% where k=3 for i,j=1,2,3, k=5 for i,j=1,2,3, k = 4 otherwise.

rho = 1025; % Water density [kg/m^3]
ULEN = 1; % WAMIT scaling factor []

T_gdf = [1 -1 -1];
Tscale = diag([T_gdf T_gdf]); % Auxiliary matrix for proprer scaling of matrices
scaleA  = [ ones(3)*3 ones(3)*4
    ones(3)*4 ones(3)*5 ];

% Extract added mass and damping
for p = 1:Nperiods
    idx = find(unique_periods == periods(p));
    Aij(dof_i(p),dof_j(p),idx) = A(p);
    Bij(dof_i(p),dof_j(p),idx) = B(p);
end

omg = (2*pi./unique_periods)'; % "omega" array, for frequency in rad/s
Nomg = length(omg); % Length of omg

% Sort with respect to frequency in increasing order
[omg_sorted, omg_idx] = sort(omg);
Aij_sorted = Aij(:,:,omg_idx);
Bij_sorted = Bij(:,:,omg_idx);
omg = omg_sorted;

% Check if the WAMIT file contains asymptotic data (i.e., coefficients for omg
% = 0 rad/s and omg = Inf).
if omg(1) == 0
    omg_asmp = 1; % Flag for indicating the presence of asymptotic data
else
    omg_asmp = 0; % No asymptotic data
end

if imemory == 0
    % If the users choses the LF + WF superposition approach, no memory effects are
    % considered and therefore there is no need for calculating the retardation
    % functions. The radiation loads considered in the equations of motions are only
    % those for omg = 0 rad/s, which may either be already provided by WAMIT or
    % extrapolated from the available data.
    if omg_asmp == 0
        A11 = zeros(6,6,1); % Matrix of added mass for ship 1
        A12 = zeros(6,6,1); % Matrix of added mass for the effects of ship 1 motions over ship 2
        A21 = zeros(6,6,1); % Matrix of added mass for the effects of ship 2 motions over ship 1
        A22 = zeros(6,6,1); % Matrix of added mass for ship 2
        for k1 = 1:6
            for k2 = 1:6
                % Extrapolate the added masses for omg = 0 rad/s
                A11(k1,k2) = interp1(omg,Aij_sorted(k1,k2,:),0,'spline','extrap');
                A12(k1,k2) = interp1(omg,Aij_sorted(k1,k2+6,:),0,'spline','extrap');
                A21(k1,k2) = interp1(omg,Aij_sorted(k1+6,k2,:),0,'spline','extrap');
                A22(k1,k2) = interp1(omg,Aij_sorted(k1+6,k2+6,:),0,'spline','extrap');
            end
        end
    elseif omg_asmp == 1
        % Take the added mass matrices for omg = 0 rad/s
        A11 = Aij_sorted(1:6,1:6,1);
        A12 = Aij_sorted(1:6,7:12,1);
        A21 = Aij_sorted(7:12,1:6,1);
        A22 = Aij_sorted(7:12,7:12,1);
    end
     
    % Scale Wamit data to SI system (Wamit axes)
    A11_dim = A11*rho .* (ULEN .^ scaleA);
    A12_dim = A12*rho .* (ULEN .^ scaleA);
    A21_dim = A21*rho .* (ULEN .^ scaleA);
    A22_dim = A22*rho .* (ULEN .^ scaleA);
    
    % Transform to Fossen axes
    A11 = Tscale*A11_dim*Tscale;
    A12 = Tscale*A12_dim*Tscale;
    A21 = Tscale*A21_dim*Tscale;
    A22 = Tscale*A22_dim*Tscale;
    
elseif imemory == 1
    % If LF and WF loads are to be considered together in the equations of
    % motions, a unified model for maneuvering and seakeeping has to be adopted
    % and therefore the radiation loads shall be represented through the convolution
    % of retardation functions. This ensures that the radiation effects due to
    % the ship motions during previous moments (i.e., the memory effects) are
    % properly considered.
    % The retardation functions are calculated from the radiation damping, as
    % described e.g. in (Journee, 1993). It is convenient for that to redefine
    % the omg array into equaly spaced intervals. Also, radiation damping for
    % high frequency values usually are not well calculated by standard hydro-
    % dynamic software, such that the tail should be approximated based on
    % available data.
    
    omg_e = linspace(omg(1),10,200); % Equaly spaced frequency, ranging from 0
    % rad/s to 10 rad/s (considered as the infinite frequency), with 200 values.
    
    delta_omg_e = omg_e(2) - omg_e(1); % Length of omg_e increment
    
    B11 = zeros(6,6,200); % Matrix of radiation damping for ship 1
    B12 = zeros(6,6,200); % Matrix of radiation damping for the effects of ship 1 motions over ship 2
    B21 = zeros(6,6,200); % Matrix of radiation damping for the effects of ship 2 motions over ship 1
    B22 = zeros(6,6,200); % Matrix of radiation damping for ship 2
    
    if omg_asmp == 0
        
        [omg_M,idx_M] = min(abs(omg_e-omg(Nomg))); % Find index of the frequency in
        % omg_e closest to omg(Nomg) (i.e., the highest available frequency
        % originally given.
        
    elseif omg_asmp == 1
        idx_M = min(abs(omg_e-3.14)); % Even when asymptotic data is available for
        % omg = Inf,it is better to considere only the data below a given limit
        % value for omg. It was decided to fiz this value at omg = 3.14 rad/s,
        % since it is not likely that WAMIT will be run for many periods between
        % 0 s and 2 s.
    end
    
    for k1 = 1:6
        for k2 = 1:6
            % Interpolate the radiation damping over the equaly spaced 
            % frequency array (until the index idx_M found above).
            for k3 = 1:idx_M
                b_aux = reshape(Bij_sorted(k1,k2,:),[1 Nomg]);
                B11(k1,k2,k3) = interp1(omg,b_aux,omg_e(k3),'spline','extrap');
                b_aux = reshape(Bij_sorted(k1,k2+6,:),[1 Nomg]);
                B12(k1,k2,k3) = interp1(omg,b_aux,omg_e(k3),'spline','extrap');
                b_aux = reshape(Bij_sorted(k1+6,k2,:),[1 Nomg]);
                B21(k1,k2,k3) = interp1(omg,b_aux,omg_e(k3),'spline','extrap');
                b_aux = reshape(Bij_sorted(k1+6,k2+6,:),[1 Nomg]);
                B22(k1,k2,k3) = interp1(omg,b_aux,omg_e(k3),'spline','extrap');
            end
            
            % According to (Journee, 1993), the tail of the radiation damping
            % curve may be approximated according to B_tail(omg) = B_max/omg^3,
            % where B_max is the damping for the highest frequency provided.
            for k3 = idx_M+1:200
                B11(k1,k2,k3) = B11(k1,k2,idx_M)/omg_e(k3)^3;
                B12(k1,k2,k3) = B12(k1,k2,idx_M)/omg_e(k3)^3;
                B21(k1,k2,k3) = B21(k1,k2,idx_M)/omg_e(k3)^3;
                B22(k1,k2,k3) = B22(k1,k2,idx_M)/omg_e(k3)^3;
            end
        end
    end
    
    % Scaling of B11, B12, B21 and B22 - adapted from (MSS, 2010)
    for k3 = 1:200        
        B11_dim = B11(:,:,k3)*rho .* omg_e(k3).* (ULEN .^ scaleA);
        B12_dim = B12(:,:,k3)*rho .* omg_e(k3).* (ULEN .^ scaleA);
        B21_dim = B21(:,:,k3)*rho .* omg_e(k3).* (ULEN .^ scaleA);
        B22_dim = B22(:,:,k3)*rho .* omg_e(k3).* (ULEN .^ scaleA);
        B11(:,:,k3) = Tscale*B11_dim*Tscale;
        B12(:,:,k3) = Tscale*B12_dim*Tscale;
        B21(:,:,k3) = Tscale*B21_dim*Tscale;
        B22(:,:,k3) = Tscale*B22_dim*Tscale;
    end   
    
    % Creation of matrices with "delta B's", that is, the difference
    % between radiation damping values for subsequent frequencies.
    delta_B11 = B11(:,:,2:200) - B11(:,:,1:199);
    delta_B12 = B12(:,:,2:200) - B12(:,:,1:199);
    delta_B21 = B21(:,:,2:200) - B21(:,:,1:199);
    delta_B22 = B22(:,:,2:200) - B22(:,:,1:199);
    
    % Initialization of matrices for storage of the retardation functions,
    % K11, K12, K21 and K22
    K11 = zeros(6,6,199);
    K12 = zeros(6,6,199);
    K21 = zeros(6,6,199);
    K22 = zeros(6,6,199);
    
    % Calculation of the retardation functions for tau = 0 (where tau is
    % the delay).
    for k1 = 1:6
        for k2 = 1:6
            b11_aux = reshape(B11(k1,k2,:),[1 200]);
            K11(k1,k2) = 2/pi*trapz(b11_aux)*delta_omg_e;
            b12_aux = reshape(B12(k1,k2,:),[1 200]);
            K12(k1,k2) = 2/pi*trapz(b12_aux)*delta_omg_e;
            b21_aux = reshape(B21(k1,k2,:),[1 200]);
            K21(k1,k2) = 2/pi*trapz(b21_aux)*delta_omg_e;
            b22_aux = reshape(B22(k1,k2,:),[1 200]);
            K22(k1,k2) = 2/pi*trapz(b22_aux)*delta_omg_e;
        end
    end
    
    % Calculation of the limit values for the integrations (Journee, 1993)
    eps = 0.01; % Factor proposed in the mentioned reference
    
    G11 = zeros(6,6);
    G12 = zeros(6,6);
    G21 = zeros(6,6);
    G22 = zeros(6,6);
    for k1 = 1:6
        for k2 = 1:6
            G11(k1,k2) = 2*sqrt((sum(abs(delta_B11(k1,k2,:))))/(pi*delta_omg_e*eps*K11(k1,k2,1)));
            G12(k1,k2) = 2*sqrt((sum(abs(delta_B12(k1,k2,:))))/(pi*delta_omg_e*eps*K12(k1,k2,1)));
            G21(k1,k2) = 2*sqrt((sum(abs(delta_B21(k1,k2,:))))/(pi*delta_omg_e*eps*K21(k1,k2,1)));
            G22(k1,k2) = 2*sqrt((sum(abs(delta_B22(k1,k2,:))))/(pi*delta_omg_e*eps*K22(k1,k2,1)));
        end
    end
    


    
%     for k1 = 1:6
%         for k2 = 1:6
%             for k3 = 1:199
%                 K11(k1,k2,k3) = K11(k1,k2,k3) + 
%             end
%         end
%     end
    
    
end








% save memory_ss A11 A12 A21 A22 B11 B12 B21 B22 omg dof_unst11 dof_unst12 dof_unst21 dof_unst22

% clear all

% References
% Journee, J. M. J. (1993) - Hydromechanics coefficients for calculating time
% domain motions of cutter suction dredges by Cummins equations. Report 968, Delft
% University of Technology, Ship Hydromechanics Laboratory. The Netherlands.
%
% MSS (2010) - Marine Systems Simulator. Viewed 24/11/2015
% http://www.marinecontrol.org.
