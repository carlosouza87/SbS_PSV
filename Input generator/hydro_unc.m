% clear all;close all;clc

rho = 1025;
ULEN = 1;

T_gdf = [1 -1 -1];
Tscale = diag([T_gdf T_gdf]);

inpt1 = ['psv.1'];
inpt2 = ['fpso.1'];

% inpt1 = [caso(k1_out).arq_isolado_vlcc '.1'];
% inpt2 = [caso(k1_out).arq_isolado_shuttle '.1'];


[periods,dof_i,dof_j,A11,B11] = textread(inpt1);
[periods,dof_i,dof_j,A22,B22] = textread(inpt2);

Nperiods = length(periods);
unique_periods = unique(periods);

% frequencies (inf is chosen as 10 rad/s)
freqs = [0 10 (2*pi./unique_periods(3:length(unique_periods)))'];

% extract added mass and damping
for p = 1:Nperiods
    idx = find(unique_periods == periods(p));
    Aij11(dof_i(p),dof_j(p),idx) = A11(p);
    Bij11(dof_i(p),dof_j(p),idx) = B11(p);
    Aij22(dof_i(p),dof_j(p),idx) = A22(p);
    Bij22(dof_i(p),dof_j(p),idx) = B22(p);
end

% sort with respect to frequency
[freqs_sorted, freq_idx] = sort(freqs);
A11 = Aij11(:,:,freq_idx);
B11 = Bij11(:,:,freq_idx);
A22 = Aij22(:,:,freq_idx);
B22 = Bij22(:,:,freq_idx);

A12 = zeros(6,6,length(freqs));
A21 = zeros(6,6,length(freqs));
B12 = zeros(6,6,length(freqs));
B21 = zeros(6,6,length(freqs));

freqs = freqs_sorted;
Nfreqs = length(freqs);

% Scaling of added mass and damping matrices (Wamit manual p. 4-3)
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

% %% Frequency-domain identification
% % system identification options
% FDIopt.OrdMax     = 20;
% FDIopt.Method     = 2;
% FDIopt.Iterations = 20;
% FDIopt.LogLin     = 1;
% FDIopt.wsFactor   = 0.1;
% FDIopt.wminFactor = 0.1;
% FDIopt.wmaxFactor = 5;
% 
% Nmax   = length(freqs);
% 
% Nf = Nmax-1;    % for WAMIT computations, remove infinite frequency data
% 
% % frequencies
% w = freqs(1:Nf)';   % does not include ininite frequency
% 
idx_M = {1,1,1,2,2,2,3,3,4,4,5,6};
idx_N = {1,3,5,2,4,6,3,5,4,6,5,6};
% 
% 
% for k1 = 1:length(idx_M),
%     
%     dof =[idx_M{k1},idx_N{k1}]   % 12 elements due to starboard-port symmetry
%     
%     % hydrodynamic raw data
%     A11_dof = reshape(A11(dof(1),dof(2),1:Nf),1,Nf)';
%     A22_dof = reshape(A22(dof(1),dof(2),1:Nf),1,Nf)';
%     B11_dof = reshape(B11(dof(1),dof(2),1:Nf),1,Nf)';
%     B22_dof = reshape(B22(dof(1),dof(2),1:Nf),1,Nf)';
%     
%     A11inf_dof = A11(dof(1),dof(2),Nmax);
%     A22inf_dof = A22(dof(1),dof(2),Nmax);
%     
%     if sum(A11(dof(1),dof(2),:)) ~= 0
%         FDIopt.AinfFlag = 1;
%         [K11rad,A11inf_hat] = FDIRadMod_alt(w,A11_dof,A11inf_dof,B11_dof,FDIopt,dof);
%     end
%     if sum(A22(dof(1),dof(2),:)) ~= 0
%         FDIopt.AinfFlag = 1;
%         [K22rad,A22inf_hat] = FDIRadMod_alt(w,A22_dof,A22inf_dof,B22_dof,FDIopt,dof);
%     end
%     
%     % compute state-space model (Ar,Br,Cr)
%     if sum(A11(dof(1),dof(2),:)) == 0
%         A11rad = 0;
%         B11rad = 0;
%         C11rad = 0;
%     else
%         [A11rad,B11rad,C11rad,D11rad] = tf2ss(K11rad.num{1},K11rad.den{1});
%     end
%    if sum(A22(dof(1),dof(2),:)) == 0
%         A22rad = 0;
%         B22rad = 0;
%         C22rad = 0;
%     else
%         [A22rad,B22rad,C22rad,D22rad] = tf2ss(K22rad.num{1},K22rad.den{1});
%    end
%     
%     A11r{dof(1),dof(2)} = A11rad;
%     B11r{dof(1),dof(2)} = B11rad;
%     C11r{dof(1),dof(2)} = C11rad;
%     D11r{dof(1),dof(2)} = 0;
%     A22r{dof(1),dof(2)} = A22rad;
%     B22r{dof(1),dof(2)} = B22rad;
%     C22r{dof(1),dof(2)} = C22rad;
%     D22r{dof(1),dof(2)} = 0;
% end

load fpsoABC

A11r = vesselABC.Ar;
B11r = vesselABC.Br;
C11r = vesselABC.Cr;
D11r = vesselABC.Dr;

A22r = vesselABC.Ar;
B22r = vesselABC.Br;
C22r = vesselABC.Cr;
D22r = vesselABC.Dr;

% symmetric elements (WAMIT and VERES)
A11r{3,1} = A11r{1,3};
B11r{3,1} = B11r{1,3};
C11r{3,1} = C11r{1,3};
D11r{3,1} = D11r{1,3};

A11r{4,2} = A11r{2,4};
B11r{4,2} = B11r{2,4};
C11r{4,2} = C11r{2,4};
D11r{4,2} = D11r{2,4};

A11r{5,1} = A11r{1,5};
B11r{5,1} = B11r{1,5};
C11r{5,1} = C11r{1,5};
D11r{5,1} = D11r{1,5};

A11r{5,3} = A11r{3,5};
B11r{5,3} = B11r{3,5};
C11r{5,3} = C11r{3,5};
D11r{5,3} = D11r{3,5};

A11r{6,2} = A11r{2,6};
B11r{6,2} = B11r{2,6};
C11r{6,2} = C11r{2,6};
D11r{6,2} = D11r{2,6};

A11r{6,4} = A11r{4,6};
B11r{6,4} = B11r{4,6};
C11r{6,4} = C11r{4,6};
D11r{6,4} = D11r{4,6};

A22r{3,1} = A22r{1,3};
B22r{3,1} = B22r{1,3};
C22r{3,1} = C22r{1,3};
D22r{3,1} = D22r{1,3};

A22r{4,2} = A22r{2,4};
B22r{4,2} = B22r{2,4};
C22r{4,2} = C22r{2,4};
D22r{4,2} = D22r{2,4};

A22r{5,1} = A22r{1,5};
B22r{5,1} = B22r{1,5};
C22r{5,1} = C22r{1,5};
D22r{5,1} = D22r{1,5};

A22r{5,3} = A22r{3,5};
B22r{5,3} = B22r{3,5};
C22r{5,3} = C22r{3,5};
D22r{5,3} = D22r{3,5};

A22r{6,2} = A22r{2,6};
B22r{6,2} = B22r{2,6};
C22r{6,2} = C22r{2,6};
D22r{6,2} = D22r{2,6};

A22r{6,4} = A22r{4,6};
B22r{6,4} = B22r{4,6};
C22r{6,4} = C22r{4,6};
D22r{6,4} = D22r{4,6};

count = 0;
dof_unst11 = [];
for k1 = 1:12
    dof =[idx_M{k1},idx_N{k1}];
    r = eig(A11r{dof(1),dof(2)});
    if sum((real(r)) > 0) > 0
        for k2 = 1:length(r)
            if real(r(k2)) > 0
                r(k2) = -real(r(k2)) + imag(r(k2))*1i;
            end
        end
        [num,den] = ss2tf(A11r{dof(1),dof(2)},B11r{dof(1),dof(2)},C11r{dof(1),dof(2)},D11r{dof(1),dof(2)});
        den = poly(r);
        [A11r{dof(1),dof(2)},B11r{dof(1),dof(2)},C11r{dof(1),dof(2)},D11r{dof(1),dof(2)}] = tf2ss(num,den);
        A11r{dof(2),dof(1)} = A11r{dof(1),dof(2)};
        B11r{dof(2),dof(1)} = B11r{dof(1),dof(2)};
        C11r{dof(2),dof(1)} = C11r{dof(1),dof(2)};
        D11r{dof(2),dof(1)} = D11r{dof(1),dof(2)};        
    end
end

count = 0;
dof_unst22 = [];
for k1 = 1:12
    dof =[idx_M{k1},idx_N{k1}];
    r = eig(A22r{dof(1),dof(2)});
    if sum((real(r)) > 0) > 0
        for k2 = 1:length(r)
            if real(r(k2)) > 0
                r(k2) = -real(r(k2)) + imag(r(k2))*1i;
            end
        end
        [num,den] = ss2tf(A22r{dof(1),dof(2)},B22r{dof(1),dof(2)},C22r{dof(1),dof(2)},D22r{dof(1),dof(2)});
        den = poly(r);
        [A22r{dof(1),dof(2)},B22r{dof(1),dof(2)},C22r{dof(1),dof(2)},D22r{dof(1),dof(2)}] = tf2ss(num,den);
        A22r{dof(2),dof(1)} = A22r{dof(1),dof(2)};
        B22r{dof(2),dof(1)} = B22r{dof(1),dof(2)};
        C22r{dof(2),dof(1)} = C22r{dof(1),dof(2)};
        D22r{dof(2),dof(1)} = D22r{dof(1),dof(2)};
    end
end

save memory_ss A11 A12 A21 A22 B11 B12 B21 B22 A11r A22r B11r B22r C11r C22r D11r D22r freqs dof_unst11 dof_unst22

% clear all
