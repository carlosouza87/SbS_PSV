function [mu] = convolution_integral(K,Tij_max,ktime,nu,t,newdt)
global variable ivalida


dt = variable.dt;
lt = variable.lt;
mu= zeros(12,1);

if ktime==1
    variable.mu = zeros (12,lt);
else
    if ivalida==0
        nu = variable.nu;
    end
    G = zeros(12,ktime);
end

if ktime > 1 && newdt == 1
    if ivalida==0
        nu(:,ktime)=variable.nu(:,ktime-1);
    end
    
    
    %Os efeitos de uma embarcação sobre a outra são considerados fazendo-se o
    %for, dentro da função, variar de k2 = 1:12
    %%%%%% Para VALIDACAO k1 e k2 foram de 1:12 para 1:1
    if newdt==1
        for k1 = 1:12
            for k2 = 1:12
                if Tij_max(k1) > t
                    for k3 = 1:ktime-1; %coloquei 2
                        %                    k3 = ktime-2;
                        %a avaliação de mu se dá até o instante k3. k3 faz o papel de
                        %"tau" do papper de Journéé e ktime indica o instante t
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


% Defina por i qual grau de liberdade deseja plotar
% i=1;
% if t>30 && newdt == 1
% figure(32)
% plot(variable.t_plot(1:ktime,1),variable.mu(i,1:ktime)','o-','LineWidth',2)
% xlabel ('tempo')
% ylabel (['mu(' num2str(i) ')'])
% grid on
% end


% i=1;
% variable.tau_plot(ktime,1)=ktime*dt;
% if newdt == 1 && ktime > lt-2
% figure(33)
% plot(variable.tau_plot(1:ktime,1),G(i,1:ktime)','o-','LineWidth',2)
% xlabel ('tau')
% ylabel (['K*nu(' num2str(i) ')'])
% grid on
% t %para mostrar tempo de simulação
% end


end

