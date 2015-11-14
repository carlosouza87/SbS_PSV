function [Xktot,Yktot,Nktot]=amarras(xg,yg,psi,Lpp,B)

global Her nlinhas

load linhadi

%xg=0;
%yg=0;
%psi=0;


%Não entendi, comentei. Se deixa ele anda em x.
% % xg=xg+0.7798154*cos(psi);
% % yg=yg+0.7798154*sin(psi);


% xa(i), ya(i) posicao absoluta da ancora
% i: amarra
% xr, yr posicao do fairlead (onde amarra a corda na embarcação) relativos ao cg


xa(1)=1269;
ya(1)=961;
xr(1)=142.1;
% yr(1)=15.4;
yr(1)=27.25;

xa(2)=1216;
ya(2)=1020;
xr(2)=139.8;
% yr(2)=16.5;
yr(2)=27.25;

xa(3)=214;
ya(3)=1490;
xr(3)=137.5;
% yr(3)=17.4;
yr(3)=27.25;

xa(4)=135;
ya(4)=1492;
xr(4)=135.2;
% yr(4)=18.4;
yr(4)=27.25;

xa(5)=56;
ya(5)=1490;
xr(5)=132.9;
% yr(5)=19.4;
yr(5)=27.25;

xa(6)=-1269;
ya(6)=961;
xr(6)=-142.1;
% yr(6)=10.1;
yr(6)=27.25;

xa(7)=-1216;
ya(7)=1020;
xr(7)=-139.8;
% yr(7)=8.4;
yr(7)=27.25;

xa(8)=-214;
ya(8)=1490;
xr(8)=-137.5;
% yr(8)=6.8;
yr(8)=27.25;

xa(9)=-135;
ya(9)=1492;
xr(9)=-135.2;
% yr(9)=4.8;
yr(9)=27.25;

xa(10)=-56;
ya(10)=1490;
xr(10)=-132.9;
% yr(10)=-4.8;
yr(10)=27.25;

xa(11)=-1269;
ya(11)=-961;
xr(11)=-142.1;
yr(11)=-27.25;

xa(12)=-1216;
ya(12)=-1020;
xr(12)=-139.8;
yr(12)=-27.25;

xa(13)=-214;
ya(13)=-1490;
xr(13)=-137.5;
yr(13)=-27.25;

xa(14)=-135;
ya(14)=-1492;
xr(14)=-135.2;
yr(14)=-27.25;

xa(15)=-56;
ya(15)=-1490;
xr(15)=-132.9;
yr(15)=-27.25;

xa(16)=1269;
ya(16)=-961;
xr(16)=142.1;
% yr(1)=15.4;
yr(16)=-27.25;

xa(17)=1216;
ya(17)=-1020;
xr(17)=139.8;
% yr(2)=16.5;
yr(17)=-27.25;

xa(18)=214;
ya(18)=-1490;
xr(18)=137.5;
% yr(3)=17.4;
yr(18)=-27.25;

xa(19)=135;
ya(19)=-1492;
xr(19)=135.2;
% yr(4)=18.4;
yr(19)=-27.25;

xa(20)=56;
ya(20)=-1490;
xr(20)=132.9;
% yr(5)=19.4;
yr(20)=-27.25;
%mantive as posições absolutas das amarras e fiz a proporção para as
%relativas.

if max(yr)> B/2 || max (xr)> Lpp/2
error ('fairlead is outside of the ship');
end

% 
% aux=0;
% for aux=6:13
% xa(aux)=0;
% ya(aux)=0;
% xr(aux)=0;
% % yr(6)=10.1;
% yr(aux)=0;
% end

nlinhas=20;

Her=zeros(1,nlinhas);
Xktot=0;
Yktot=0;
Xkrtot=0;
Ykrtot=0;
Nktot=0;
for k=1:nlinhas
% xf(k) posicao absoluta do fairlead
 xf(k)=xg + xr(k)*cos(psi)-yr(k)*sin(psi);%lembrar que o navio gira em torno do CG, e essa linha é de vital importancia para nmodificarmos o lr a medida que o
%navio se move.
 yf(k)=yg + xr(k)*sin(psi)+yr(k)*cos(psi);
 lr(k)=sqrt((xf(k)-xa(k))^2+(yf(k)-ya(k))^2);

 %He e l representam uma relação entre a força exercida pela corda ao ser
 %tracionada, amarrada na embarcação e em um outro ponto.
 
 rompeu=0;
 if lr(k)>max(l)
  rompeu=k
 end

 if rompeu==k
     Her(k)=max(He);
 else  
  i=1;
  while lr(k)>l(i)  
   i=i+1;
  end
 
  if i==1
   Her(k)=0; %ou seja, no caso do lr ser menor que já o primeiro i
  else
   Her(k)=He(i-1)+(lr(k)-l(i-1))*(He(i)-He(i-1))/(l(i)-l(i-1));
%interpolação por reta entre os pontos (He(i-1);l(i-1)) e (He(i);l(i)).
  end
 end
 
 

 Xk(k)=Her(k)*(xa(k)-xf(k))/lr(k); %ele projeta a força Her(k) em x, com base em lrx(k)/lr(k).
 Yk(k)=Her(k)*(ya(k)-yf(k))/lr(k); 
 Xkr(k)=Xk(k)*cos(psi)+Yk(k)*sin(psi);%ele passa a força em x da corda  para o referencial da embarcação.
 Ykr(k)=Yk(k)*cos(psi)-Xk(k)*sin(psi);
 Nk(k)=-Xkr(k)*yr(k)+Ykr(k)*xr(k); %Nk(k) é o torque no eixo Z a linha de amarração k;
 Xktot=Xktot+Xk(k); %acha a força resultante no referencial absoluto.
 Yktot=Yktot+Yk(k);
 Xkrtot=Xkrtot+Xkr(k);%acha a força resultante no referencial da embarcação.
 Ykrtot=Ykrtot+Ykr(k);
 Nktot=Nktot+Nk(k);
end



