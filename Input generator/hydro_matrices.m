clear all;close all;clc

% It is admitted that the WAMIT files brings either both asymptotic values (omg
% = 0 rad/s and omg = Inf) or none of them.


caseid = 'conjunto';
inpt1 = [caseid '.1']; % Input file

% Define the values for the dimensionalization
rho = 1025;
ULEN = 1;

T_gdf = [1 -1 -1];
Tscale = diag([T_gdf T_gdf]); % Matrix for proprer scaling of values

% Read the input file, assigning each column to the corresponding array
[periods,dof_i,dof_j,A,B] = textread(inpt1);

Nperiods = length(periods); % Number of periods, considering repeated
unique_periods = unique(periods); % Array of non-repeating periods

% Extract added mass and damping
for p = 1:Nperiods
    idx = find(unique_periods == periods(p));
    Aij(dof_i(p),dof_j(p),idx) = A(p);
    Bij(dof_i(p),dof_j(p),idx) = B(p);
end

omg = (2*pi./unique_periods)'; % "omega" array, for frequency in rad/s
Nomg = length(omg); % Current number of elements of omg

% Sort with respect to frequency in increasing order
[omg_sorted, omg_idx] = sort(omg);
Aij_sorted = Aij(:,:,omg_idx);
Bij_sorted = Bij(:,:,omg_idx);
omg = omg_sorted;

% A treatment is now performed for ensuring that:
% - the first values of added mass and radiation damping correspond to omg = 
%   0 rad/s;
% - the radiation damping is truncated to zero after the last given value.
%   For that, the highest value of omg is incremented and appended to omg,
%   such that the corresponding damping is set to zero and kept constant
%   thenceforth. The added mass for the incremented omg is the same as for
%   the previous value, and is also kept constant thenceforth.
% - omg = 10 rad/s is assigned as the assymptotic infinite frequency, and
%   is therefore the last value of the omg array.
% 
% The array of frequencies is therefore updated for proper truncation of 
% added mass and radiation damping. This means:
% - appending a omg = 0 rad/s at the beginning of the array, if the original
%   data does not contain it;
% - appending an incremental value after the highest available frequency
%   (disconsidering the infinite frequency, if available;
% - appending omg = 10 rad/s at the end of the array, corresponding to
%   infinite frequency (or setting omg = inf to omg = 10 rad/s, if omg =
%   inf is available in the original data).
 
d_omg = 0.01; % Incremental frequency (rad/s)

% If there is no value assigned for omg = 0 rad/s in the  original data,
% append omg = 0 rad/s as the first element.
if omg_sorted(1) == 0
    omg_asmp = 1; % Flag for indicating that the WAMIT file contains asymptotic 
    % data (i.e., omg = 0 rad/s and omg = Inf)  
    omg = [omg(1:Nomg-1) 10]; % In this case, a new frequency, the only 
    % change is the replacement of the last value (originally omg = Inf) by 
    % omg = 10 rad/s.                        
else
    omg_asmp = 0; % No asymptotic data in WAMIT file
    omg = [0 omg omg(Nomg)+d_omg 10]; % In this case, omg = 0 rad/s is appended
    % to the beginning of the frequency array. Also, the last frequency is
    % incremented with d_omg and appended to the end of the array. Finally, omg
    % = 10 rad/s is appended to the end of the array.
end
Nomg = length(omg); % Update Nomg

% Now the matrices with added mass and radiation damping coefficientes are
% updated according to the new omg arrays

% Declare new matrices for storing the values
A11 = zeros(6,6,Nomg);
A12 = zeros(6,6,Nomg);
A21 = zeros(6,6,Nomg);
A22 = zeros(6,6,Nomg);
B11 = zeros(6,6,Nomg);
B12 = zeros(6,6,Nomg);
B21 = zeros(6,6,Nomg);
B22 = zeros(6,6,Nomg);

% 
if omg_asmp == 0
    % First elements (omg = 0 rad/s) for added masses are the same as the first
    % in the original files. For radiation damping, they are forced to be zero.
    A11(:,:,1) = Aij_sorted(1:6,1:6,1);
    A12(:,:,1) = Aij_sorted(1:6,7:12,1);
    A21(:,:,1) = Aij_sorted(7:12,1:6,1);
    A22(:,:,1) = Aij_sorted(7:12,7:12,1);
    B11(:,:,1) = zeros(6,6);
    B12(:,:,1) = zeros(6,6);
    B21(:,:,1) = zeros(6,6);
    B22(:,:,1) = zeros(6,6);
    % The next elements from the original data are now appended to both the added
    % mass and radiation damping matrices. Since the original omg array was
    % increased in 3 elements, the range for the new A and B matrices is 2:Nomg-2.     
    A11(:,:,2:Nomg-2) = Aij_sorted(1:6,1:6,:);
    A12(:,:,2:Nomg-2) = Aij_sorted(1:6,7:12,:);
    A21(:,:,2:Nomg-2) = Aij_sorted(7:12,1:6,:);
    A22(:,:,2:Nomg-2) = Aij_sorted(7:12,7:12,:);
    B11(:,:,2:Nomg-2) = Bij_sorted(1:6,1:6,:);
    B12(:,:,2:Nomg-2) = Bij_sorted(1:6,7:12,:);
    B21(:,:,2:Nomg-2) = Bij_sorted(7:12,1:6,:);
    B22(:,:,2:Nomg-2) = Bij_sorted(7:12,7:12,:);
    % For the added mass, the value corresponding to the incremented frequency
    % equals the previous one (i.e., the last value in the original file). For
    % the radiation damping, it is already assigned to zero.
    A11(:,:,Nomg-1) = Aij_sorted(1:6,1:6,Nomg-3);
    A12(:,:,Nomg-1) = Aij_sorted(1:6,7:12,Nomg-3);
    A21(:,:,Nomg-1) = Aij_sorted(7:12,1:6,Nomg-3);
    A22(:,:,Nomg-1) = Aij_sorted(7:12,7:12,Nomg-3);
    B11(:,:,Nomg-1) = zeros(6,6);
    B12(:,:,Nomg-1) = zeros(6,6);
    B21(:,:,Nomg-1) = zeros(6,6);
    B22(:,:,Nomg-1) = zeros(6,6);
    % Finally, the last value for the added mass is again set equal the
    % previous one. For the radiation damping, it is again set as zero.
    A11(:,:,Nomg) = Aij_sorted(1:6,1:6,Nomg-3);
    A12(:,:,Nomg) = Aij_sorted(1:6,7:12,Nomg-3);
    A21(:,:,Nomg) = Aij_sorted(7:12,1:6,Nomg-3);
    A22(:,:,Nomg) = Aij_sorted(7:12,7:12,Nomg-3);
    B11(:,:,Nomg) = zeros(6,6);
    B12(:,:,Nomg) = zeros(6,6);
    B21(:,:,Nomg) = zeros(6,6);
    B22(:,:,Nomg) = zeros(6,6); 
else
    % In the case the WAMIT data contains values for asymptoticaly zero and 
    % infinite frequencies, the only worry is to ensure that the the 
    % radiation damping is zero for omg = 0 rad/s and for omg = Inf. All
    % the other elements (including added mass) are kept as provided by 
    % WAMIT.
    A11(:,:,1:Nomg) = Aij_sorted(1:6,1:6,1:Nomg);
    A12(:,:,1:Nomg) = Aij_sorted(1:6,7:12,1:Nomg);
    A21(:,:,1:Nomg) = Aij_sorted(7:12,1:6,1:Nomg);
    A22(:,:,1:Nomg) = Aij_sorted(7:12,7:12,1:Nomg);
    B11(:,:,1) = zeros(6,6);
    B12(:,:,1) = zeros(6,6);
    B21(:,:,1) = zeros(6,6);
    B22(:,:,1) = zeros(6,6);
    B11(:,:,Nomg) = zeros(6,6);
    B12(:,:,Nomg) = zeros(6,6);
    B21(:,:,Nomg) = zeros(6,6);
    B22(:,:,Nomg) = zeros(6,6);
end

% Scaling of added mass and damping matrices (p. 4-3 of Wamit manual v. 6.2s)
% Aij = Aij' * rho * ULEN^k
% Bij = Bij' * rho * w * ULEN^k
% where k=3 for i,j=1,2,3, k=5 for i,j=1,2,3, k = 4 otherwise.
scaleA  = [ ones(3)*3 ones(3)*4
    ones(3)*4 ones(3)*5 ];

for w = 1:Nomg
    % scale Wamit data to SI system (Wamit axes)
    A11_dim = A11(:,:,w)*rho .* (ULEN .^ scaleA);
    A12_dim = A12(:,:,w)*rho .* (ULEN .^ scaleA);
    A21_dim = A21(:,:,w)*rho .* (ULEN .^ scaleA);
    A22_dim = A22(:,:,w)*rho .* (ULEN .^ scaleA);
    B11_dim = B11(:,:,w)*rho .* omg(w) ...
        .* (ULEN .^ scaleA);
    B12_dim = B12(:,:,w)*rho .* omg(w) ...
        .* (ULEN .^ scaleA);
    B21_dim = B21(:,:,w)*rho .* omg(w) ...
        .* (ULEN .^ scaleA);
    B22_dim = B22(:,:,w)*rho .* omg(w) ...
        .* (ULEN .^ scaleA);
    % transform to Fossen axes
    A11(:,:,w) = Tscale*A11_dim*Tscale;
    A12(:,:,w) = Tscale*A12_dim*Tscale;
    A21(:,:,w) = Tscale*A21_dim*Tscale;
    A22(:,:,w) = Tscale*A22_dim*Tscale;
    B11(:,:,w) = Tscale*B11_dim*Tscale;
    B12(:,:,w) = Tscale*B12_dim*Tscale;
    B21(:,:,w) = Tscale*B21_dim*Tscale;
    B22(:,:,w) = Tscale*B22_dim*Tscale;
end

for k1 = 1:Nomg
    a11(k1) = A22(1,1,k1);
    a22(k1) = A22(2,2,k1);
    a33(k1) = A22(3,3,k1);
    a44(k1) = A22(4,4,k1);
    a55(k1) = A22(5,5,k1);
    a66(k1) = A22(6,6,k1);    
    b11(k1) = B22(1,1,k1);
    b22(k1) = B22(2,2,k1);
    b33(k1) = B22(3,3,k1);
    b44(k1) = B22(4,4,k1);
    b55(k1) = B22(5,5,k1);
    b66(k1) = B22(6,6,k1); 
end

figure(1)
plot(omg,a11)
grid on
figure(2)
plot(omg,a22)
grid on
figure(3)
plot(omg,a33)
grid on
figure(4)
plot(omg,a44)
grid on
figure(5)
plot(omg,a55)
grid on
figure(6)
plot(omg,a66)
grid on
figure(7)
plot(omg,b11)
grid on
figure(8)
plot(omg,b22)
grid on
figure(9)
plot(omg,b33)
grid on
figure(10)
plot(omg,b44)
grid on
figure(11)
plot(omg,b55)
grid on
figure(12)
plot(omg,b66)
grid on


% Colocar um if aqui, e ativar apenas se o usu�rio quiser. Incluir fun��o
% do andrey




% save memory_ss A11 A12 A21 A22 B11 B12 B21 B22 omg dof_unst11 dof_unst12 dof_unst21 dof_unst22

% clear all
