function [t,y,nn]=hrkdif(y,dt,t,nn)
nn=t;
newdt=1;
[dy,y,nsys,nn]=eqsim(y,t,newdt,nn); %1
% nn
newdt=0; %Só considera o primeiro eqsim
h=dt/2;
% if t>50
% y
% end

for k1=1:nsys
	sy(k1)=y(k1);
	y0(k1)=dy(k1);
	y(k1)=h*dy(k1)+y(k1);
end

t=t+h;

[dy,y,nsys,nn]=eqsim(y,t,newdt,nn); %2
% nn
for k1=1:nsys
	y1(k1)=dy(k1);
	y(k1)=sy(k1)+h*dy(k1);
end

[dy,y,nsys,nn]=eqsim(y,t,newdt,nn); %3
% nn
for k1=1:nsys
	y2(k1)=dy(k1);
	y(k1)=sy(k1)+dt*dy(k1);
end;

t=t+h;

[dy,y,nsys,nn]=eqsim(y,t,newdt,nn); %4
% nn
h=h/3; %esse h/3 que corresponde a dt/6 é usado somente no loop abaixo

for k1=1:nsys
	prt1=2*(y1(k1)+y2(k1));
	prt2=y0(k1)+dy(k1);
	y(k1)=sy(k1)+h*prt1+h*prt2;
end