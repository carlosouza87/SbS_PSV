function [K,dK,Tij,Nt,B,A] = retardation_function(B,A)
global data variable

%recebe-se um A e B que são saidas do wamit e o devolvem interpolados para
%frequencias linearmente espaçadas
%Pegou-se agora um intervalo de freqs que tivesse, desconsiderando o Ainf,
%a mesma abrangência de frequências mas que seja linearmente espaçado.

lt = variable.lt;
dt = variable.dt;
freqs = data.hydro.freqs;
N = length(freqs);
w =linspace(freqs(1),freqs(N),N);% arbitrei prolongar o range de frequências em duas vezes
dw= w(2)-w(1);
w_aux = linspace(freqs(N)+dw,2*freqs(N)+dw,N); %intervalo até Ainfinito
w = [w w_aux];
w = w';
freqs = freqs';

B_freqs=zeros(length(freqs),1);
B_esp=zeros(6,6,length(w));

A_freqs=zeros(length(freqs),1);
A_esp=zeros(6,6,length(freqs));

dB=zeros(6,6,length(w));

Tij=zeros(6);
Nt=zeros(6);

%Agora interpola-se para cada grau de liberdade expressos por k1 e k2
%qual o valor nas novas frequencias
for k1 = 1:6
    for k2 = 1:6
        %B_freqs representa um vetor auxiliar com os
        %valores do elemento B(k1,k2) em todas as frequências.
        for k3 = 1:N; 
            B_freqs(k3) = B(k1,k2,k3);
            A_freqs(k3) = A(k1,k2,k3);
         end
        
        %B_esp é o vetor de B linearmente espaçado, resultado da
        %interpolação de w(k3), ou seja, todas as novas frequências.
        for k3 = 1:N;
        B_esp(k1,k2,k3) = interp1(freqs, B_freqs,w(k3),'spline');
        A_esp(k1,k2,k3) = interp1(freqs, A_freqs,w(k3),'spline');
        
        end
        for k3 = 1:N; %tail of curve
          B_esp(k1,k2,N+k3) = B_esp(k1,k2,N)/w_aux(k3)^3; % a cauda da curva será "completada" com o ultimo valor fornecido pelo WAMIT dividido por w^3, tal que esse w é a posição no eixo x.
        end
    end
end


B=B_esp;
A=A_esp;


% i=2;
% A_plot=zeros(size(A,3),1);
% figure(5)
% A_plot(:,1)= A(i,i,:);
% plot(w(1:N),A_plot,'o-','LineWidth',2)
% xlabel ('frequencia w')
% ylabel (['A(' num2str(i) ',' num2str(i) ')' ]);
% grid on
% 
% B_plot=zeros(size(B,3),1);
% figure(6)
% B_plot(:,1)= B(i,i,:);
% plot(w,B_plot,'o-','LineWidth',2)
% xlabel ('frequencia w')
% ylabel (['B(' num2str(i) ',' num2str(i) ')']);
% grid on


%criação do vetor delta B, que representa a diferença de amortecimento
%entre duas frequencias consecutivas, para o instante inicial

for k1 = 1:6
    for k2 = 1:6
        for k3 = 1:length(w);
            if k3==1
            dB(k1,k2,k3) = B(k1,k2,k3);    
            else
            %é importante guardar dessa forma,
            %pois posteriormente, para o calculo do (tempo limite) Tij deve-se saber o grau i e j (k1,k2)
            %e também entre quais duas frequências é feita a diferença
            dB(k1,k2,k3) = B(k1,k2,k3)-B(k1,k2,k3-1); 
            end
            
        end
    end
 end
K=zeros(6,6,lt);
 %Cálculo do K(k1,k2,1) correspondente ao instante inicial
 for k1 = 1:6
    for k2 = 1:6   
     K_aux=0;
        for k3 = 1:length(w) %seria o omega do infinite added mass
     K_aux = K_aux + B(k1,k2,k3)*dw;            
        end
     K(k1,k2,1)=2/pi*K_aux;    
    end
 end
 
dB_soma = zeros(6,6); 

for k1 = 1:6
    for k2 = 1:6
        for k3 = 1:length(w);
        dB_soma(k1,k2) = dB_soma(k1,k2) + abs(dB(k1,k2,k3));
        end
        e=0.010; %segundo papper de Journéé
    Tij(k1,k2)=2*( dB_soma(k1,k2) /(pi*dw*e*K(k1,k2,1)) )^0.5;
    Tij(k1,k2)=fix(Tij(k1,k2));
    end
end
Tij=abs(Tij); %o que fazer quando Tij é complexo, causado por (Kij(k1,k2,0)<0)^0.5, resultado da soma dos B(k1,k2,k3)*dw que 
%às vezes possui termos negativos, quando um grau de liberdade 'joga' energia sobre o outro?
% Nt=abs(Nt);
%Cálculo da função de retardo (retardation functions) até o instante Tij,sendo o time step igual a dt

% Nt_biggest = max(max(Nt));
% aux = K;
% K=zeros(6,6,Nt_biggest);
% K(:,:,1)= aux;


 for k1 = 1:6
    for k2 = 1:6
        
      %K_aux representa a somatória de n=1 a Nw como está no Papper de Journéé
        
        for tau_calc= 2:lt %É necessário calcular para um grau (k1,k2) até o tempo total, exceto no instante inicial que já foi calculado.
            tau = tau_calc*dt;
            K_aux=0;
            for k3 = 1:length(w);
                if k3==1
            K_aux = K_aux + dB(k1,k2,k3)/dw*(cos(w(k3)*tau));
            
                else
            K_aux = K_aux + dB(k1,k2,k3)/dw*(cos(w(k3)*tau)-cos(w(k3-1)*tau ));
                end
            end
       
        K(k1,k2,tau_calc) = 2/(pi*tau^2)*K_aux + 2/(pi*tau)*B(k1,k2,length(w))*sin(w(k3)*tau);

        end
    end
 end

%criação do vetor delta K, que representa a diferença da retardation
%function entre dois instantes de tempo consecutivos, para cada grau (k1,k2) de liberdade
dK= zeros(6,6,size(K,3));


for k1=1:6
        for k2=1:6
            for k3=1:size(K,3)
               if k3==1
                dK(k1,k2,1)= K(k1,k2,k3);
               else
                 dK(k1,k2,k3)= K(k1,k2,k3)-K(k1,k2,k3-1);
               end
           end
        end
end


data.hydro.w = w;

% i=5;
% figure(33)
% t_plot= (1:size(K,3))*dt;
% t_plot=t_plot';
% K_plot(:,1)= K(i,i,:);
% plot(t_plot,K_plot(1:length(t_plot),1),'o-','LineWidth',2)
% xlabel ('tempo')
% ylabel (['K(' num2str(i) ',' num2str(i) ')'])
% grid on

end

