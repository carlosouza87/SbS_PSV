function [l1,l2,l3,l4,l5] = shipdraw(Xg,Yg,psi,L,B,redfac,frmori,colmat)
% [L1,L2,L3,L4,L5] = shipdraw(Xg,Yg,psi,L,B,redfac,colmat)
% Desenha um navio sobre o gráfico atual. Xg e Yg são as coordenadas 
% globais do centro de gravidade  do navio. Psi é o aproamento (rad), L o
% comprimento, B a boca, 0 < redfac < 1 é a escala de redução do desenho, 
% frmori é a orientação do sistema de coordenadas (1 NWU, -1 para NED)* e 
% colmat é uma matriz 1X3 com  as  componentes RGB para a cor das linhas. 
% Se esse parâmetro não for dado, as linhas serão pretas por padrão (i.e., 
% colmat = [0 0 0]).
% * NWU: North-West-UP. NED: North-East-Down
% 
% Carlos Souza, 25/01/2011 - Universidade de São Paulo

if frmori == -1
    coordned = [0 1 0;1 0 0;0 0 -1]*[Xg;Yg;psi];
    Xg = coordned(1);
    Yg = coordned(2);
    psi = coordned(3);
    psi = psi + pi/2;
end
if nargin == 6
    colmat = [0 0 0];
end
Pl(1,:) = redfac*[-L/2,frmori*-B/2];
Pl(2,:) = redfac*[0.7*L/2,frmori*-B/2];
Pl(3,:) = redfac*[L/2,frmori*0];
Pl(4,:) = redfac*[0.7*L/2,frmori*B/2];
Pl(5,:) = redfac*[-L/2,frmori*B/2];

Jpsi = [cos(psi) -sin(psi);sin(psi) cos(psi)];
for k = 1:5
    P(:,k) = [Xg;Yg] + Jpsi*Pl(k,:)';    
end
P = P';
l1 = line([P(1,1),P(2,1)],[P(1,2),P(2,2)],'Color',colmat);
l2 = line([P(2,1),P(3,1)],[P(2,2),P(3,2)],'Color',colmat);
l3 = line([P(3,1),P(4,1)],[P(3,2),P(4,2)],'Color',colmat);
l4 = line([P(4,1),P(5,1)],[P(4,2),P(5,2)],'Color',colmat);
l5 = line([P(5,1),P(1,1)],[P(5,2),P(1,2)],'Color',colmat);