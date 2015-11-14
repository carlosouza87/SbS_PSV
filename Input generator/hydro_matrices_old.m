% clear all;close all;clc
% clear variables;close all;clc

rho = 1025;
ULEN = 1;

T_gdf = [1 -1 -1];
Tscale = diag([T_gdf T_gdf]);

% fileinput  = input('File name: ','s');
% if strcmp(fileinput,'')
%     vesselfile = filename;
% else
%     vesselfile = fileinput;       % default value
% end

%inpt1 = [file_hydro_matrices '.1'];
inpt1 = ['conjunto.1'];
[periods,dof_i,dof_j,A,B] = textread(inpt1);

Nperiods = length(periods);
unique_periods = unique(periods);

% frequencies (inf is set as 10 rad/s)
freqs = [0 10 (2*pi./unique_periods(3:length(unique_periods)))']; % CORRIGIR - a(w=10) = a(w_max), e b(w=10) = 0

% extract added mass and damping
for p = 1:Nperiods
    idx = find(unique_periods == periods(p));
    Aij(dof_i(p),dof_j(p),idx) = A(p);
    Bij(dof_i(p),dof_j(p),idx) = B(p);
end

% sort with respect to frequency
[freqs_sorted, freq_idx] = sort(freqs);
Aij_sorted = Aij(:,:,freq_idx);
Bij_sorted = Bij(:,:,freq_idx);

A11 = Aij_sorted(1:6,1:6,:);
A12 = Aij_sorted(1:6,7:12,:);
A21 = Aij_sorted(7:12,1:6,:);
A22 = Aij_sorted(7:12,7:12,:);

B11 = Bij_sorted(1:6,1:6,:);
B12 = Bij_sorted(1:6,7:12,:);
B21 = Bij_sorted(7:12,1:6,:);
B22 = Bij_sorted(7:12,7:12,:);

freqs = freqs_sorted;
Nfreqs = length(freqs);

% Scaling of added mass and damping matrices (Wamit manual p. 4-3) CHECAR p. p/ v 6.2
% Aij = Aij' * rho * ULEN^k
% Bij = Bij' * rho * w * ULEN^k
% where k=3 for i,j=1,2,3, k=5 for i,j=1,2,3, k = 4 otherwise.
scaleA  = [ ones(3)*3 ones(3)*4
    ones(3)*4 ones(3)*5 ];

for w = 1:Nfreqs,
    % scale Wamit data to SI system (Wamit axes)
    A11_dim = A11(:,:,w)*rho .* (ULEN .^ scaleA);
    A12_dim = A12(:,:,w)*rho .* (ULEN .^ scaleA);
    A21_dim = A21(:,:,w)*rho .* (ULEN .^ scaleA);
    A22_dim = A22(:,:,w)*rho .* (ULEN .^ scaleA);
    B11_dim = B11(:,:,w)*rho .* freqs(w) ...
        .* (ULEN .^ scaleA);
    B12_dim = B12(:,:,w)*rho .* freqs(w) ...
        .* (ULEN .^ scaleA);
    B21_dim = B21(:,:,w)*rho .* freqs(w) ...
        .* (ULEN .^ scaleA);
    B22_dim = B22(:,:,w)*rho .* freqs(w) ...
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

% Colocar um if aqui, e ativar apenas se o usuário quiser. Incluir função
% do andrey




save memory_ss A11 A12 A21 A22 B11 B12 B21 B22 A11r A12r A21r A22r B11r B12r B21r B22r C11r C12r C21r C22r D11r D12r D21r D22r freqs dof_unst11 dof_unst12 dof_unst21 dof_unst22

% clear all
