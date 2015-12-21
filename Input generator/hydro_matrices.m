% Reads added mass and radiation damping from WAMIT .1 file. Then, output retardation
% functions (and the respective time limits) and infinite-frequency added mass
% matrices, if Cummins equations are considered, or zero-frequency added mass
% matrices, if the LF+WF motions superposition approach is adopted. For either
% case, matrices have indeces 11 (corresponding to ship 1), 12 (ship 1 over ship
% 2), 21 (ship 2 over ship 1) and 22 (ship 2).

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

T_gdf = [1 -1 -1];
Tscale = diag([T_gdf T_gdf]); % Matrix for transforming data into Fossen axes
scaleA  = [ ones(3)*3 ones(3)*4
    ones(3)*4 ones(3)*5 ]; % Matrix for scaling WAMIT data

% Extract added mass and damping
for p = 1:Nperiods
    idx = find(unique_periods == periods(p));
    Aij(dof_i(p),dof_j(p),idx) = A(p);
    Bij(dof_i(p),dof_j(p),idx) = B(p);
end

omg = (2*pi./unique_periods)'; % "omega" array, for frequency in rad/s
Nomg = length(omg); % Number of frequencies considered

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

if isimtype == 1
    % If LF and WF loads are to be considered together in the equations of
    % motions, a unified model for maneuvering and seakeeping must be adopted
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
        
        B11(:,:,k3) = B11_dim;
        B12(:,:,k3) = B12_dim;
        B21(:,:,k3) = B21_dim;
        B22(:,:,k3) = B22_dim;
        
        %        % Transform to Fossen axes
        %         B11(:,:,k3) = Tscale*B11*Tscale;
        %         B12(:,:,k3) = Tscale*B12*Tscale;
        %         B21(:,:,k3) = Tscale*B21*Tscale;
        %         B22(:,:,k3) = Tscale*B22*Tscale;
    end
    
    % Creation of matrices with "delta B's", that is, the difference
    % between radiation damping values for subsequent frequencies.
    delta_B11 = B11(:,:,2:200) - B11(:,:,1:199);
    delta_B12 = B12(:,:,2:200) - B12(:,:,1:199);
    delta_B21 = B21(:,:,2:200) - B21(:,:,1:199);
    delta_B22 = B22(:,:,2:200) - B22(:,:,1:199);
    
    % Initialization of matrices for storage of the retardation functions,
    % K11_0, K12_0, K21_0 and K22_0
    K11_0 = zeros(6,6,1);
    K12_0 = zeros(6,6,1);
    K21_0 = zeros(6,6,1);
    K22_0 = zeros(6,6,1);
    
    % Calculation of the retardation functions for tau = 0 s (where tau is
    % the delay).
    for k1 = 1:6
        for k2 = 1:6
            b11_aux = reshape(B11(k1,k2,:),[1 200]);
            K11_0(k1,k2) = 2/pi*trapz(b11_aux)*delta_omg_e;
            b12_aux = reshape(B12(k1,k2,:),[1 200]);
            K12_0(k1,k2) = 2/pi*trapz(b12_aux)*delta_omg_e;
            b21_aux = reshape(B21(k1,k2,:),[1 200]);
            K21_0(k1,k2) = 2/pi*trapz(b21_aux)*delta_omg_e;
            b22_aux = reshape(B22(k1,k2,:),[1 200]);
            K22_0(k1,k2) = 2/pi*trapz(b22_aux)*delta_omg_e;
        end
    end
    
    %    % Calculation of the limit values for the integrations (Journee, 1993)
    %    eps = 0.01; % Factor proposed in the mentioned reference
    %
    %    T11 = zeros(6,6);
    %    T12 = zeros(6,6);
    %    T21 = zeros(6,6);
    %    T22 = zeros(6,6);
    %    for k1 = 1:6
    %        for k2 = 1:6
    %            T11(k1,k2) = 2*sqrt((sum(abs(delta_B11(k1,k2,:))))/(pi*delta_omg_e*eps*K11_0(k1,k2,1)));
    %            T12(k1,k2) = 2*sqrt((sum(abs(delta_B12(k1,k2,:))))/(pi*delta_omg_e*eps*K12_0(k1,k2,1)));
    %            T21(k1,k2) = 2*sqrt((sum(abs(delta_B21(k1,k2,:))))/(pi*delta_omg_e*eps*K21_0(k1,k2,1)));
    %            T22(k1,k2) = 2*sqrt((sum(abs(delta_B22(k1,k2,:))))/(pi*delta_omg_e*eps*K22_0(k1,k2,1)));
    %        end
    %    end
    % Fix all limit values to 120 s
    T11 = ones(6,6)*120;
    T12 = ones(6,6)*120;
    T21 = ones(6,6)*120;
    T22 = ones(6,6)*120;
    
    % Calculation of retardation functions for tau > 0 s (Journee, 1993, sec. 4.4)
    for k1 = 1:6
        for k2 = 1:6
            tau11 = dt:dt:T11(k1,k2);
            K11(k1,k2).tau = [0 tau11];
            K_aux = zeros(1,length(tau11));
            for k4 = 1:199
                K_aux = K_aux + (2./(pi*tau11.^2)).*(delta_B11(k1,k2,k4)/delta_omg_e*(cos(omg_e(k4+1)*tau11)-cos(omg_e(k4)*tau11)));
            end
            m1 = 2./(pi*tau11);
            m2 = B11(k1,k2,200)*sin(omg_e(200)*tau11);
            K_aux = K_aux + m1.*m2;
            K11(k1,k2).K = [K11_0(k1,k2) K_aux];
            K11(k1,k2).T = T11(k1,k2);
            
            tau12 = dt:dt:T12(k1,k2);
            K12(k1,k2).tau = [0 tau12];
            K_aux = zeros(1,length(tau12));
            for k4 = 1:199
                K_aux = K_aux + (2./(pi*tau12.^2)).*(delta_B12(k1,k2,k4)/delta_omg_e*(cos(omg_e(k4+1)*tau12)-cos(omg_e(k4)*tau12)));
            end
            m1 = 2./(pi*tau12);
            m2 = B12(k1,k2,200)*sin(omg_e(200)*tau12);
            K_aux = K_aux + m1.*m2;
            K12(k1,k2).K = [K12_0(k1,k2) K_aux];
            K12(k1,k2).T = T12(k1,k2);
            
            tau21 = dt:dt:T21(k1,k2);
            K21(k1,k2).tau = [0 tau21];
            K_aux = zeros(1,length(tau21));
            for k4 = 1:199
                K_aux = K_aux + (2./(pi*tau21.^2)).*(delta_B21(k1,k2,k4)/delta_omg_e*(cos(omg_e(k4+1)*tau21)-cos(omg_e(k4)*tau21)));
            end
            m1 = 2./(pi*tau21);
            m2 = B21(k1,k2,200)*sin(omg_e(200)*tau21);
            K_aux = K_aux + m1.*m2;
            K21(k1,k2).K = [K21_0(k1,k2) K_aux];
            K21(k1,k2).T = T21(k1,k2);
            
            tau22 = dt:dt:T22(k1,k2);
            K22(k1,k2).tau = [0 tau22];
            K_aux = zeros(1,length(tau22));
            for k4 = 1:199
                K_aux = K_aux + (2./(pi*tau22.^2)).*(delta_B22(k1,k2,k4)/delta_omg_e*(cos(omg_e(k4+1)*tau22)-cos(omg_e(k4)*tau22)));
            end
            m1 = 2./(pi*tau22);
            m2 = B22(k1,k2,200)*sin(omg_e(200)*tau22);
            K_aux = K_aux + m1.*m2;
            K22(k1,k2).K = [K22_0(k1,k2) K_aux];
            K22(k1,k2).T = T22(k1,k2);
        end
    end
    
    % Calculation of infinite added mass matrices, if not provided by WAMIT
    if omg_asmp == 0
        % Coefficients for the first frequency provided in the WAMIT file
        % will be used for the calculation.
        omg_fix = omg(10);
        A11 = Aij_sorted(1:6,1:6,10);
        A12 = Aij_sorted(1:6,7:12,10);
        A21 = Aij_sorted(7:12,1:6,10);
        A22 = Aij_sorted(7:12,7:12,10);
        % Scale Wamit data to SI system (Wamit axes)
        A11 = A11*rho .* (ULEN .^ scaleA);
        A12 = A12*rho .* (ULEN .^ scaleA);
        A21 = A21*rho .* (ULEN .^ scaleA);
        A22 = A22*rho .* (ULEN .^ scaleA);
        
        %         % Transform to Fossen axes
        %         A11 = Tscale*A11*Tscale;
        %         A12 = Tscale*A12*Tscale;
        %         A21 = Tscale*A21*Tscale;
        %         A22 = Tscale*A22*Tscale;
        
        A11_inf = zeros(6,6);
        A12_inf = zeros(6,6);
        A21_inf = zeros(6,6);
        A22_inf = zeros(6,6);
        
        OMG = omg_e(200);
        for k1 = 1:6
            for k2 = 1:6
                tau11 = K11(k1,k2).tau;
                k11 = K11(k1,k2).K; % Lower case k11 for not overwriting K11
                lk = length(k11);
                Int11 = 0;
                for k3 = 2:200
                    dK = k11(k3) - k11(k3-1);
                    Int11 = Int11 + 1/OMG^2*dK/dt*(sin(OMG*k3*dt)-sin(OMG*(k3-1)*dt));
                end
                Int11 = Int11 + 1/OMG*(k11(1)-k11(lk)*cos(OMG*200*dt));
                A11_inf(k1,k2) = A11(k1,k2) + 1/omg_fix*Int11;
                
                tau12 = K12(k1,k2).tau;
                k12 = K12(k1,k2).K; % Lower case k12 for not overwriting K12
                lk = length(k12);
                Int12 = 0;
                for k3 = 2:200
                    dK = k12(k3) - k12(k3-1);
                    Int12 = Int12 + 1/OMG^2*dK/dt*(sin(OMG*k3*dt)-sin(OMG*(k3-1)*dt));
                end
                Int12 = Int12 + 1/OMG*(k12(1)-k12(lk)*cos(OMG*200*dt));
                A12_inf(k1,k2) = A12(k1,k2) + 1/omg_fix*Int12;
                
                tau21 = K21(k1,k2).tau;
                k21 = K21(k1,k2).K; % Lower case k21 for not overwriting K21
                lk = length(k21);
                Int21 = 0;
                for k3 = 2:200
                    dK = k21(k3) - k21(k3-1);
                    Int21 = Int21 + 1/OMG^2*dK/dt*(sin(OMG*k3*dt)-sin(OMG*(k3-1)*dt));
                end
                Int21 = Int21 + 1/OMG*(k21(1)-k21(lk)*cos(OMG*200*dt));
                A21_inf(k1,k2) = A21(k1,k2) + 1/omg_fix*Int21;
                
                tau22 = K22(k1,k2).tau;
                k22 = K22(k1,k2).K; % Lower case k22 for not overwriting K22
                lk = length(k22);
                Int22 = 0;
                for k3 = 2:200
                    dK = k22(k3) - k22(k3-1);
                    Int22 = Int22 + 1/OMG^2*dK/dt*(sin(OMG*k3*dt)-sin(OMG*(k3-1)*dt));
                end
                Int22 = Int22 + 1/OMG*(k22(1)-k22(lk)*cos(OMG*200*dt));
                A22_inf(k1,k2) = A22(k1,k2) + 1/omg_fix*Int22;
            end
        end
        
    end
elseif isimtype == 2
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
                A_interp = reshape(Aij_sorted(k1,k2,:),Nomg,1);
                A11(k1,k2) = interp1(omg,A_interp,0,'spline','extrap');
                A_interp = reshape(Aij_sorted(k1,k2+6,:),Nomg,1);
                A12(k1,k2) = interp1(omg,A_interp,0,'spline','extrap');
                A_interp = reshape(Aij_sorted(k1+6,k2,:),Nomg,1);
                A21(k1,k2) = interp1(omg,A_interp,0,'spline','extrap');
                A_interp = reshape(Aij_sorted(k1+6,k2+6,:),Nomg,1);
                A22(k1,k2) = interp1(omg,A_interp,0,'spline','extrap');
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
    A11 = A11*rho .* (ULEN .^ scaleA);
    A12 = A12*rho .* (ULEN .^ scaleA);
    A21 = A21*rho .* (ULEN .^ scaleA);
    A22 = A22*rho .* (ULEN .^ scaleA);
    
    %     % Transform to Fossen axes
    %     A11 = Tscale*A11*Tscale;
    %     A12 = Tscale*A12*Tscale;
    %     A21 = Tscale*A21*Tscale;
    %     A22 = Tscale*A22*Tscale;
end

% Save relevant variables to be used in the simulations
if isimtype == 1
    save('hydro_data','A11_inf','A12_inf','A21_inf','A22_inf','K11','K12','K21','K22','dt')
elseif isimtype == 2
    save('hydro_data','A11','A12','A21','A22')
end

% References
% Journee, J. M. J. (1993) - Hydromechanics coefficients for calculating time
% domain motions of cutter suction dredges by Cummins equations. Report 968, Delft
% University of Technology, Ship Hydromechanics Laboratory. The Netherlands.
%
% MSS (2010) - Marine Systems Simulator. Viewed 24/11/2015
% http://www.marinecontrol.org.
