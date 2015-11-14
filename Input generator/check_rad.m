% check_rad - plots radiation coefficients read from WAMIT file
% works for WAMIT files corresponding to either an isolate ship or to a set
% of two ships. In the latter case, the coefficients for each ship and for
% the influence of a ship over the other are plotted.
%
% The routine must be run in the same folder where the WAMIT files are placed.
%
% Last edited on 11/11/2015 by Carlos Souza - University of Sao Paulo


clear all;close all;clc

% caseid = input('Input case ID: ');
caseid = 'conjunto';

rho = 1025;
ULEN = 1;

T_gdf = [1 -1 -1];
Tscale = diag([T_gdf T_gdf]);

inpt1 = [caseid '.1'];
[periods,dof_i,dof_j,A,B] = textread(inpt1);

% Determine number of ships in the file. If there is one ship only, the number
% of DOFs is 6. For two ships, the DOFs of the second ship range from 7 to 12.
if max(dof_i) == 6
    nships = 1;
else
    nships = 2;
end

Nperiods = length(periods); % Number of periods, considering repeated
unique_periods = unique(periods); % List with unique (non-repeating) periods

% Extract added mass and damping
for p = 1:Nperiods
    idx = find(unique_periods == periods(p));
    Aij(dof_i(p),dof_j(p),idx) = A(p);
    Bij(dof_i(p),dof_j(p),idx) = B(p);
end

omg = (2*pi./unique_periods)'; % "Omega" array, for frequency in rad/s
Nomg = length(omg); % Number of frequencies in the array

% Sort with respect to frequency in increasing order
[omg_sorted, omg_idx] = sort(omg);
Aij_sorted = Aij(:,:,omg_idx);
Bij_sorted = Bij(:,:,omg_idx);

%% Ploting data
% dof_inp = input('Input a pair of DOF inside brackets (e.g. [1 3]): ');

%dof_inp = input('Input a pair of DOF (e.g. [1 3]), or ''all'' for ploting all combinations: ');
%if dof_inp == 'all'
%    sym_xy = input('Symmetry in XY plane? (y/n)');
%    sym_yz = input('Symmetry in YZ plane? (y/n)');
%    sym_xz = input('Symmetry in XZ plane? (y/n)');
%end

DOFs = [1 1;1 3;1 5;2 2;2 4;2 6;3 1;3 3;3 5;4 2;4 4;4 6;5 1;5 3;5 5;6 2;6 4;6 6];
% DOFs = [1 2;1 4;1 6;2 1;2 3;2 5;3 2;3 4;3 6];
for k_dof = 1:size(DOFs,1)
    close all   
    dof_inp = DOFs(k_dof,:);
    if nships == 1
        a = zeros(1,Nomg);
        b = zeros(1,Nomg);
        for k_omg = 1:Nomg
            a(k_omg) = Aij_sorted(dof_inp(1),dof_inp(2),k_omg); % Added mass for the DOF of interest
            b(k_omg) = Bij_sorted(dof_inp(1),dof_inp(2),k_omg); % Radiation damping for the DOF of interest
        end
        figure(1)
        plot(omg,a)
        title_str = ['Added mass - ' caseid ' - a' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Added mass (non-dimensional)')
        grid on
        figure(2)
        plot(omg,b)
        title_str = ['Radiation damping - ' caseid ' - b' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Radiation damping (non-dimensional)')
        grid on
    else
        % When two ships are present, plots are generated for each ship (11 and 22)
        % and for their mutual influence (12) and (21).
        a11 = zeros(1,Nomg); % Added mass, ship 1
        b11 = zeros(1,Nomg); % Radiation damping, ship 1
        a12 = zeros(1,Nomg); % Added mass, ship 1 over ship 2
        b12 = zeros(1,Nomg); % Radiation damping, ship 1 over ship 2
        a21 = zeros(1,Nomg); % Added mass, ship 2 over ship 1
        b21 = zeros(1,Nomg); % Radiation damping, ship 2 over ship 1
        a22 = zeros(1,Nomg); % Added mass, ship 2
        b22 = zeros(1,Nomg); % Radiation damping, ship 2
        for k_omg = 1:Nomg
            a11(k_omg) = Aij_sorted(dof_inp(1),dof_inp(2),k_omg); % Added mass for the DOF of interest, ship 1
            b11(k_omg) = Bij_sorted(dof_inp(1),dof_inp(2),k_omg); % Radiation damping for the DOF of interest, ship 1
            a12(k_omg) = Aij_sorted(dof_inp(1),dof_inp(2)+6,k_omg); % Added mass for the DOF of interest, ship 1 over ship 2
            b12(k_omg) = Bij_sorted(dof_inp(1),dof_inp(2)+6,k_omg); % Radiation damping for the DOF of interest, ship 1 over ship 2
            a21(k_omg) = Aij_sorted(dof_inp(1)+6,dof_inp(2),k_omg); % Added mass for the DOF of interest, ship 2 over ship 1
            b21(k_omg) = Bij_sorted(dof_inp(1)+6,dof_inp(2),k_omg); % Radiation damping for the DOF of interest, ship 2 over ship 1
            a22(k_omg) = Aij_sorted(dof_inp(1)+6,dof_inp(2)+6,k_omg); % Added mass for the DOF of interest, ship 2
            b22(k_omg) = Bij_sorted(dof_inp(1)+6,dof_inp(2)+6,k_omg); % Radiation damping for the DOF of interest, ship 2
        end
        figure(1)
        subplot(2,2,1)
        plot(omg,a11)
        title_str = ['Added mass - ' caseid ' - ship 1 - a' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Added mass (non-dimensional)')
        grid on
        subplot(2,2,2)
        plot(omg,a12)
        title_str = ['Added mass - ' caseid ' - ship 1 over ship 2 - a' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Added mass (non-dimensional)')
        grid on
        subplot(2,2,3)
        plot(omg,a21)
        title_str = ['Added mass - ' caseid ' - ship 2 over ship 1 - a' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Added mass (non-dimensional)')
        grid on
        subplot(2,2,4)
        plot(omg,a22)
        title_str = ['Added mass - ' caseid ' - ship 2 - a' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Added mass (non-dimensional)')
        grid on
        figure(2)
        subplot(2,2,1)
        plot(omg,b11)
        title_str = ['Radiation damping - ' caseid ' - ship 1 - b' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Radiation damping (non-dimensional)')
        grid on
        subplot(2,2,2)
        plot(omg,b12)
        title_str = ['Radiation damping - ' caseid ' - ship 1 over ship 2 - b' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Radiation damping (non-dimensional)')
        grid on
        subplot(2,2,3)
        plot(omg,b21)
        title_str = ['Radiation damping - ' caseid ' - ship 2 over ship 1 - b' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Radiation damping (non-dimensional)')
        grid on
        subplot(2,2,4)
        plot(omg,b22)
        title_str = ['Radiation damping - ' caseid ' - ship 2 - b' num2str(dof_inp(1)) num2str(dof_inp(2))];
        title(title_str)
        xlabel('\omega (rad/s)' )
        ylabel('Radiation damping (non-dimensional)')
        grid on
    end
    pause
end