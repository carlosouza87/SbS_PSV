function [A_inf] = infinite_added_mass(A,K,dK)
global data variable
lt = variable.lt;
freqs = data.hydro.freqs;
dt = variable.dt;
w = data.hydro.w;
omega = w(length(w)); %infinite frequency
A_inf = zeros(6);

for k1=1:6
    for k2=1:6
    K_aux=zeros(6);
        for k3=1:lt
        K_aux(k1,k2)= K_aux(k1,k2)+ dK(k1,k2,k3)/dt*(sin(omega*k3*(k3*dt))-sin(omega*(k3-1)*dt));
        end
        K_aux(k1,k2)= 1/omega^2*K_aux(k1,k2) + 1/omega*( K(k1,k2,1)- K(k1,k2,lt)*cos(omega*lt*dt));
%          K_aux(k1,k2)= 1/omega^2*K_aux(k1,k2) + 1/omega*( K(k1,k2,1)- K(k1,k2,Nt(k1,k2))*cos(omega*Nt(k1,k2)*dt)); 
        A_inf(k1,k2)= A(k1,k2,length(freqs))+1/omega*K_aux(k1,k2);
    end
end

end