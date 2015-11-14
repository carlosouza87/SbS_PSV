function [tau_m1,tau_m2,sgm,swl_exc,mbl_exc] = mooring(eta1,eta2,lines,stiff,swl,mbl,newdt) %problema
% function [tau_m1,tau_m2,sgm,swl_exc,mbl_exc] = mooring(eta1,eta2,pos1,pos2,lines,stiff,swl,mbl)
% Calculates mooring loads and verifies wether SWL and MBL are exceeded or
% not.
%
% Inputs:
% eta1 - position and Euler vectors of the ship 1
% eta2 - position and Euler vectors of the ship 2
% lines - lines geometrical properties
% e.g.: lines = [lengthrope1 lengthrope2 ... lengthropeN
%                   Drope1     Drope2    ...   Drope3
%                   X1rope1    X1rope2   ...   X1ropeN
%                   Y1rope1    Y1rope2   ...   Y1ropeN
%                   Z1rope1    Z1rope2   ...   Z1ropeN
%                   X2rope1    X2rope2   ...   X2ropeN
%                   Y2rope1    Y2rope2   ...   Y2ropeN
%                   Z2rope1    Z2rope2   ...   Z2ropeN ]
%
% where D is the diameter of the rope and (X1, Y1 and Z1) and (X2, Y2 and
% Z2) are the connection points in the  ships 1 and 2, respectivelly.
% Units are in S.I.
% stiff - lines stiffnesses %matriz repetida do programa do carlos eduardo
% e.g.: stiff = [0        0.02    ...     0.10
%                T1(0)  T1(0.02)  ...    T1(0.10)
%                T2(0)  T2(0.02)  ...    T2(0.10)
%                  .      .        .       .
%                  .      .        .       .
%                  .      .        .       .
%                TN(0)  TN(0.02)  ...    TN(0.10)]
% where the first row contains rope elongation, and TX is the rope
% tension?) acredito o correto ser: rope FORCE em newtons.
% for each of the respective elongations of the rope X.
% swl - Service Working Load (N)
% mbl - Minimum Breaking Load (N)
% 
% Outputs:
% tau_m1,tau_m2 - 6X1 vectors with the loads induced by the mooring lines
% on each ship [N, Nm]
% sgm - vector with the tensions 
% swl_exc,mbl_exc - flags to check whether swl and mbl are exceeded (1) or
% not (0)
% 
% Carlos Souza, 24/01/2011 - Universidade de São Paulo

nlin = size(lines,2);%acho que o dois da função size se referencia a colunas
tau_m1 = zeros(6,nlin);
tau_m2 = zeros(6,nlin);
sgm = zeros(nlin,1); %vetor com as tensoes;
swl_exc = 0;
mbl_exc = 0;

for k1 = 1:nlin
    phi1 = eta1(4);
    theta1 = eta1(5);
    psi1 = eta1(6);
%     psi1=0.26;
    
    
    Rth1 = [cos(psi1)*cos(theta1) -sin(psi1)*cos(phi1)+cos(psi1)*sin(theta1)*sin(phi1) sin(psi1)*sin(phi1)+cos(psi1)*cos(phi1)*sin(theta1);
        sin(psi1)*cos(theta1) cos(psi1)*cos(phi1)+sin(phi1)*sin(theta1)*sin(psi1) -cos(psi1)*sin(phi1)+sin(theta1)*sin(psi1)*cos(phi1);
        -sin(theta1) cos(theta1)*sin(phi1) cos(theta1)*cos(phi1)];
    
    %     Tth1 = [1 sin(phi1)*tan(theta1) cos(phi1)*tan(theta1);
    %         0 cos(phi1) -sin(phi1);
    %         0 sin(phi1)/cos(theta1) cos(phi1)/cos(theta1)];
    
    %     Jth1 = [Rth1 zeros(3);zeros(3) Tth1];
    
    %     P1 = eta1 + Jth1*[lines(3,k1);lines(4,k1);lines(5,k1);0;0;0];
    P1 = eta1 + [Rth1*[lines(3,k1);lines(4,k1);lines(5,k1)]; 0; 0; 0];
    X1 = P1(1);
    Y1 = P1(2);
    Z1 = P1(3);
    
    Pproa1 = eta1 + [Rth1*[166.5;27.25;-6.3]; 0; 0; 0]; %no caso  geral, seria loa/2;
    
    
    Ppopa1 = eta1 + [Rth1*[-166.5;27.25;-6.3]; 0; 0; 0]; %em geral, seria -loa/2;  
%     Ppopa1 = eta1 + [Rth1*[-166.5;0;-6.3]; 0; 0; 0]  
    
    %phi2 = 0.26;
    phi2 = eta2(4);
    theta2 = eta2(5);
    psi2 = eta2(6);
%     psi2 = 0.26;
%     
    Rth2 = [cos(psi2)*cos(theta2) -sin(psi2)*cos(phi2)+cos(psi2)*sin(theta2)*sin(phi2) sin(psi2)*sin(phi2)+cos(psi2)*cos(phi2)*sin(theta2);
        sin(psi2)*cos(theta2) cos(psi2)*cos(phi2)+sin(phi2)*sin(theta2)*sin(psi2) -cos(psi2)*sin(phi2)+sin(theta2)*sin(psi2)*cos(phi2);
        -sin(theta2) cos(theta2)*sin(phi2) cos(theta2)*cos(phi2)];
    
    %     Tth2 = [1 sin(phi2)*tan(theta2) cos(phi2)*tan(theta2);
    %         0 cos(phi2) -sin(phi2);
    %         0 sin(phi2)/cos(theta2) cos(phi2)/cos(theta2)];
    
    %     Jth2 = [Rth2 zeros(3);zeros(3) Tth2];
    
    %     P2 = eta2 + Jth2*[lines(6,k1);lines(7,k1);lines(8,k1);0;0;0];
    P2 = eta2 + [Rth2*[lines(6,k1);lines(7,k1);lines(8,k1)]; 0; 0; 0];
    X2 = P2(1);
    Y2 = P2(2);
    Z2 = P2(3);
    
    
    Pproa2 = eta2 + [Rth2*[44;-9.5;-1.4]; 0; 0; 0]; %no caso do PSV, em geral, seria loa/2;
    Yproa2 = Pproa2(2);

    
    Ppopa2 = eta2 + [Rth2*[-44;-9.5;-1.4]; 0; 0; 0]; %em geral, seria -loa/2;  
    Ypopa2 = Ppopa2(2);
   
    %proa
    dir1_proa= (Ppopa1-Pproa1)/norm(Ppopa1-Pproa1); %versor com a direção da embarcação 1, e sentido apontando para Ppopa1.
    delta_proa = Pproa2 - Pproa1; %vetor de diferença entre a proa das duas embarcações, apontando para a Pproa2.
    Yproa1_comp = Pproa1(2) + norm(delta_proa)* dir1_proa(2); 
    
    %popa
    dir1_popa= (Pproa1-Ppopa1)/norm(Pproa1-Ppopa1); %versor com a direção da embarcação 1, e sentido apontando para Pproa1.
    delta_popa = Ppopa2 - Ppopa1; %vetor de diferença entre a popa das duas embarcações, apontando para a Ppopa2.
    Ypopa1_comp = Ppopa1(2) + norm(delta_popa)* dir1_popa(2); 

    if( Yproa2 - Yproa1_comp <= 0)
        error ('choque entre as embarcações');
    end
    
    if( Ypopa2 - Ypopa1_comp <= 0) %a expressão a direita da subtração representa uma correção para transportarmos o ponto Ypopa1 
%que está no extremo até a mesma cota em X do ponto de comparação Ypopa2.
        error ('choque entre as embarcações');
    end
    
    d = sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);% distancia em aprox. 6.2 corda rompe
%     if(k1==1)
%     d
%      end
%     if (newdt==1)
%         xxxx=1;
%     end
    direc1 = [(X1-X2)/d;(Y1-Y2)/d;(Z1-Z2)/d];%criou um vetor com as mesma direções de d mas cujo modulo é 1;.
    if d > lines(1,k1)
        eps = (d-lines(1,k1))/lines(1,k1);
        %         T = interp1(stiff(1,:),stiff(k1+1,:),eps,'spline');
        T = interp1(stiff(1,:),stiff(k1+1,:),eps,'linear','extrap');
        if T < 0 %afinal nao é uma mola
            T = 0;
        end
        sgm(k1) = T/(lines(2,k1)^2*pi/4); %acho que não é utilizado no programa
        %         sgm(k1) = lines(2,k1)*eps;
        %         T = sgm(k1)*lines(3,k1)^2*pi/4;
    else
        sgm(k1) = 0;
        T = 0;
    end
    % Rope force projections (ship 1)
    Fx1 = -direc1(1)*T;
    Fy1 = -direc1(2)*T;
    Fz1 = -direc1(3)*T;
    % Force in body fixed frame (ship 1)
    F1 = Rth1\[Fx1;Fy1;Fz1];
    
    % Moment(ship 1)
    M1(1,1) = (lines(4,k1)*F1(3)-lines(5,k1)*F1(2));%faz produto vetorial para tirar momento mas eixos são ortogonais sin sempre 1.
    M1(2,1) = (lines(5,k1)*F1(1)-lines(3,k1)*F1(3));
    M1(3,1) = (lines(3,k1)*F1(2)-lines(4,k1)*F1(1));
    tau_m1(:,k1) = [F1;M1];
    % Rope force projections (ship 2)
    Fx2 = direc1(1)*T;
    Fy2 = direc1(2)*T;
    Fz2 = direc1(3)*T;
    % Force in body fixed frame (ship 2)
    F2 = Rth2\[Fx2;Fy2;Fz2];
    
    % Moment (ship 2)
    M2(1,1) = (lines(7,k1)*F2(3)-lines(8,k1)*F2(2));
    M2(2,1) = (lines(8,k1)*F2(1)-lines(6,k1)*F2(3));
    M2(3,1) = (lines(6,k1)*F2(2)-lines(7,k1)*F2(1));
    tau_m2(:,k1) = [F2;M2]; %um é força outro momento.
    if mbl_exc ~= 1
        mbl_exc = T > mbl; %aqui compara se rompeu ou não.
    end
    if swl_exc ~= 1
        swl_exc = T > swl;
    end
 %eles param de aparecer quando mbl_exc é um vetor de 4 linhas iguais a 1.   
%  mbl_exc
%  swl_exc


% if( eta2(4)> 0)
% eta1
% end

% if( mbl_exc==1 )
% error ('estourou');
% end
  if( eta2(4)> 0.52) %30 graus
error ('emborcando');
% mbl_exc ;t; 
  end
end


