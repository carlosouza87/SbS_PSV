function [mu] = convolution_integral(K_tot,nu,t,ktime)
% function [mu] = convolution_integral(K_tot,nu,t,ktime,newdt)
% global variable
% 
% dt = variable.dt;
dt = 0.1;
% lt = variable.lt;

% Vector for memory effects
mu = zeros(12,1);

% if ktime==1
%     variable.mu = zeros (12,lt);
% else
%     nu = variable.nu;
%     G = zeros(12,ktime);
% end

% if ktime > 1 && newdt == 1
if ktime > 1
    for k1 = 1:12
        for k2 = 1:12
            K = K_tot(k1,k2).K;
            T = K_tot(k1,k2).T;
            % While the current simulation time (t) is shorter than the
            % duration (t) of the retardation function (K), the integrand  
            % length is limited to the current simulation time index (ktime), 
            % such that its size matches the current size of the velocities
            % vector (nu). When the simulation time surpasses T the
            % integrand length takes length T/dt + 1 and the retardation
            % function is therefore evaluated for the whole duration of K.
            if t < T
                Int = zeros(1,ktime);
                for k3 = 1:ktime
                    Int(k3) = K(k3)*nu(k2,ktime-k3+1); % Integrand
                end
                mu(k1) = trapz(Int)*dt; % Convolution integral
            else
                k_fin = T/dt + 1; % End index
                Int = zeros(1,k_fin); 
                for k3 = 1:k_fin 
                   Int(k3) = K(k3)*nu(k2,k_fin-k3+1); % Integrand
                end
                mu(k1) = trapz(Int)*dt; % Convolution integral
            end
        end
    end
    
    %     variable.mu(:,ktime) =  mu(:,1);
    %  variable.mu(:,ktime) = variable.mu(:,ktime-1)+ mu(:,1);
end




