function [mu] = convolution_integral(K,tau_max,nu,t,ktime,newdt)
global variable

dt = variable.dt;
lt = variable.lt;

% Vector for memory effects
mu = zeros(12,1);

if ktime==1
    variable.mu = zeros (12,lt);
else
    nu = variable.nu;
    G = zeros(12,ktime);
end

if ktime > 1 && newdt == 1
    if ivalida==0
        nu(:,ktime)=variable.nu(:,ktime-1);
    end    
    
    if newdt==1
        for k1 = 1:12
            for k2 = 1:12
                if Tij_max(k1) > t
                    for k3 = 1:ktime-1; %coloquei 2
                        %                    k3 = ktime-2;
                        %a avaliacao de mu se da ate o instante k3. k3 faz o papel de
                        %"tau" do paper de Journee e ktime indica o instante t
                        %     mu(k1) = mu(k1) + ( (K(k1,k2,k3)*etap(k2,ktime-k3))+(K(k1,k2,k3+1)*etap(k2,ktime-k3-1)) )*dt/2; %posso considerar dtau=dt?
                        mu(k1) = mu(k1) +  K(k1,k2,k3)*nu(k2,ktime-k3)*dt;
                        G(k1,k3) =  K(k1,k2,k3)*nu(k2,ktime-k3);
                    end
                else
                    for k3 = 1 : Tij_max(k1)/dt
                        mu(k1) = mu(k1) +  K(k1,k2,k3)*nu(k2,ktime-k3)*dt;
                        G(k1,k3) =  K(k1,k2,k3)*nu(k2,ktime-k3);
                    end
                end
            end
        end
    end
    
    variable.mu(:,ktime) =  mu(:,1);
    
    %  variable.mu(:,ktime) = variable.mu(:,ktime-1)+ mu(:,1);
end
variable.t_plot(ktime,1)=ktime*dt;


end

