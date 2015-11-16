% clear all;close all;clc
% clear variables;close all;clc


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
if omg_sorted(1) ~= 0
    omg_zero = 0; % Flag for indicating that omg = 0 rad/s was not found in 
                  % the original data
    omg = [0 omg_sorted];  
else
    omg_zero = 1; % omg = 0 rad/s was found in the original data
end
Nomg = length(omg); % Update Nomg

% If the file brings data for infinite frequency (T = 0 s), it will be assigned 
% as "Inf" in the array of frequencies. It should then be replaced by 10 rad/s
% Therefore, the last element is replaced, the incremented frequency is
% appended and finally 10 is appended.
if omg_inf == 1
    omg = [omg(1:Nomg-1) omg(Nomg-1)+d_omg 10];                                                 
else
% Otherwise, the last element of omg is incremented and appended to the 
% array, and 10 is then appended at the end.
    omg = [omg omg(Nomg)+d_omg 10];
end
    
Nomg = length(omg); % Update Nomg

A11 = zeros(6,6,Nomg);
A12 = zeros(6,6,Nomg);
A21 = zeros(6,6,Nomg);
A22 = zeros(6,6,Nomg);
B11 = zeros(6,6,Nomg);
B12 = zeros(6,6,Nomg);
B21 = zeros(6,6,Nomg);
B22 = zeros(6,6,Nomg);




if omg_zero == 0
    A11(:,:,1) = Aij_sorted(1:6,1:6,1);
    A12(:,:,1) = Aij_sorted(1:6,7:12,1);
    A21(:,:,1) = Aij_sorted(7:12,1:6,1);
    A22(:,:,1) = Aij_sorted(7:12,7:12,1);
    B11(:,:,1) = zeros(6,6,Nomg);
    B12(:,:,1) = zeros(6,6,Nomg);
    B21(:,:,1) = zeros(6,6,Nomg);
    B22(:,:,1) = zeros(6,6,Nomg);
    if omg_inf == 0
        A11(:,:,2:Nomg-2) = Aij_sorted(1:6,1:6,:);
        A12(:,:,2:Nomg-2) = Aij_sorted(1:6,7:12,:);
        A21(:,:,2:Nomg-2) = Aij_sorted(7:12,1:6,:);
        A22(:,:,2:Nomg-2) = Aij_sorted(7:12,7:12,:);    
        B11(:,:,2:Nomg-2) = Bij_sorted(1:6,1:6,:);
        B12(:,:,2:Nomg-2) = Bij_sorted(1:6,7:12,:);
        B21(:,:,2:Nomg-2) = Bij_sorted(7:12,1:6,:);
        B22(:,:,2:Nomg-2) = Bij_sorted(7:12,7:12,:); 
        A11(:,:,Nomg-1) = Aij_sorted(1:6,1:6,Nomg-2);
        A12(:,:,Nomg-1) = Aij_sorted(1:6,7:12,Nomg-2);
        A21(:,:,Nomg-1) = Aij_sorted(7:12,1:6,Nomg-2);
        A22(:,:,Nomg-1) = Aij_sorted(7:12,7:12,Nomg-2);    
        B11(:,:,Nomg-1) = zeros(6,6,Nomg);
        B12(:,:,Nomg-1) = zeros(6,6,Nomg);
        B21(:,:,Nomg-1) = zeros(6,6,Nomg);
        B22(:,:,Nomg-1) = zeros(6,6,Nomg);
        A11(:,:,Nomg) = Aij_sorted(1:6,1:6,Nomg-2);
        A12(:,:,Nomg) = Aij_sorted(1:6,7:12,Nomg-2);
        A21(:,:,Nomg) = Aij_sorted(7:12,1:6,Nomg-2);
        A22(:,:,Nomg) = Aij_sorted(7:12,7:12,Nomg-2);    
        B11(:,:,Nomg) = zeros(6,6,Nomg);
        B12(:,:,Nomg) = zeros(6,6,Nomg);
        B21(:,:,Nomg) = zeros(6,6,Nomg);
        B22(:,:,Nomg) = zeros(6,6,Nomg);
    else
        A11(:,:,2:Nomg-2) = Aij_sorted(1:6,1:6,:);
        A12(:,:,2:Nomg-2) = Aij_sorted(1:6,7:12,:);
        A21(:,:,2:Nomg-2) = Aij_sorted(7:12,1:6,:);
        A22(:,:,2:Nomg-2) = Aij_sorted(7:12,7:12,:);  
    end
else
end

A11(:,:,2:Nomg-1) = Aij_sorted(1:6,1:6,:);
A12(:,:,2:Nomg-1) = Aij_sorted(1:6,7:12,:);
A21(:,:,2:Nomg-1) = Aij_sorted(7:12,1:6,:);
A22(:,:,2:Nomg-1) = Aij_sorted(7:12,7:12,:);



B11(:,:,2:Nomg-1) = Bij_sorted(1:6,1:6,:);

B12(:,:,2:Nomg-1) = Bij_sorted(1:6,7:12,:);

B21(:,:,2:Nomg-1)= Bij_sorted(7:12,1:6,:);

B22(:,:,2:Nomg-1)= Bij_sorted(7:12,7:12,:);

% Define the added mass for omg = 0 rad/s equal to the first value provided by WAMIT
A11(:,:,1) = A11(:,:,2);
A12(:,:,1) = A12(:,:,2);
A21(:,:,1) = A21(:,:,2);
A22(:,:,1) = A22(:,:,2);

% Define the added mass for omg = 10 rad/s equal to the last value (highest frequency)
A11(:,:,Nomg) = A11(:,:,Nomg-1);
A12(:,:,Nomg) = A12(:,:,Nomg-1);
A21(:,:,Nomg) = A21(:,:,Nomg-1);
A22(:,:,Nomg) = A22(:,:,Nomg-1);

% Define the radiation damping for omg = 0 rad/s and for omg = 10 rad/s equal to zero
B11(:,:,1) = 0;
B12(:,:,1) = 0;
B21(:,:,1) = 0;
B22(:,:,1) = 0;
B11(:,:,Nomg) = 0;
B12(:,:,Nomg) = 0;
B21(:,:,Nomg) = 0;
B22(:,:,Nomg) = 0;


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

% Colocar um if aqui, e ativar apenas se o usu�rio quiser. Incluir fun��o
% do andrey




% save memory_ss A11 A12 A21 A22 B11 B12 B21 B22 omg dof_unst11 dof_unst12 dof_unst21 dof_unst22

% clear all
