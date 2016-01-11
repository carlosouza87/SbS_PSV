function [mu] = convolution_integral(K_tot,nu,t,ktime)
% function [mu] = convolution_integral(K_tot,nu,t,ktime)
% Calculate the convolution term in Cummins equation, given the retardation
% function K_tot, the current velocity vector nu, the time instant t and
% the time step ktime.
% 
% Input:
% K_tot - previously generated structure with retardation functions and 
% their respective duration.
% nu - [2*ndof x l] velocity vector for both ships. The number of columns,
% l, equals:
%  - ktime, if t is lower than the retardation functions length, ltau
%  - ltau, otherwise
% It should be noted that if the ship keeps a mean velocity of motion in
% some DOF for a long time, this mean velocity should be subtracted from nu
% before using it in the function.
% t - current time instant.
% ktime - current time step.
% 
% Output:
% mu - [ndof x 1] vector with the convolution integral calculated for the
% current time step.

dt = 0.1;
mu = zeros(12,1);
% if ktime > 1 && newdt == 1
T = K_tot(1,1).T;
if ktime > 1
    if t < T
        nu_m = mean(nu,2);
        delta_nu = nu - repmat(nu_m,1,ktime);
    else
        k_fin = T/dt; % End index
        Nw = 61;
        wts = [1/(2*Nw) repmat(1/Nw,1,Nw) 1/(2*Nw)];
        nu_m = zeros(12,k_fin+1);
        for k3 = 1:12
            nu_m(k3,:) = mean_nu(nu(k3,ktime-k_fin:ktime),wts);
        end
        delta_nu = nu(:,ktime-k_fin:ktime) - nu_m;
    end    
    
    for k1 = 1:12
        for k2 = 1:12
            K = K_tot(k1,k2).K;            
            % While the current simulation time (t) is shorter than the
            % duration (t) of the retardation function (K), the integrand
            % length is limited to the current simulation time index (ktime),
            % such that its size matches the current size of the velocities
            % vector (nu). When the simulation time surpasses T the
            % integrand length takes length T/dt + 1 and the retardation
            % function is therefore evaluated for the whole duration of K.
            if t < T
                Int = K(1:ktime).*fliplr(delta_nu(k2,:));
                mu(k1) = mu(k1) + trapz(Int)*dt; % Convolution integral
            else
                Int = K.*fliplr(delta_nu(k2,1:k_fin+1));
                mu(k1) = mu(k1) + trapz(Int)*dt; % Convolution integral
            end
        end
    end
    
end