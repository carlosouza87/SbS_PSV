function [etap,nup] = eqmot_memory(eta,nu,rg1,rg2,Mrb,Ainf,G,mu,tau,newdt)
% function [etap,nup] = eqmot_memory(eta,nu,rg1,rg2,Mrb,Ainf,G,mu,tau)
% Performs calculation of the equations of motion, with consideration of
% fluid memory effects. Inputs are:
%
% eta - vector of positions (surge, sway) and heading, earth-fixed reference system
% nu - vector of velocities (surge, sway and yaw), ship-fixed reference system
% rg - position of center of gravity
% Mrb - Rigid body inertia matrix
% Ainf - Infinity frequency added inertia matrix
% G - restoring matrix
% mu - memory terms vector (convolution integral)
% tau - vector with external forces and moment acting on the ship [N,N,Nm]
%
% Outputs are:
%
% etap - vector with time derivatives of eta
% nup - vector with time derivatives of nu
%
% Carlos Souza, 26/11/2011 - Universidade de São Paulo

%% Nested rigid-body centripetal and Coriolis matrix generator function
    function Crb = genCrb(Mrb,nu,rg)
        xg = rg(1);
        yg = rg(2);
        zg = rg(3);
        
        m = Mrb(1,1);
        Ixz = -Mrb(4,6);
        Ixy = -Mrb(4,5);
        Iyz = -Mrb(5,6);
        Ix = Mrb(4,4);
        Iy = Mrb(5,5);
        Iz = Mrb(6,6);
        
        u = nu(1);
        v = nu(2);
        w = nu(3);
        p = nu(4);
        q = nu(5);
        r = nu(6);
        
        Crb = [0 0 0 m*(yg*q+zg*r) -m*(xg*q-w) -m*(xg*r+v);
            0 0 0 -m*(yg*p+w) m*(zg*r+xg*p) -m*(yg*r-u);
            0 0 0 -m*(zg*p-v) -m*(zg*q+u) m*(xg*p+yg*q);
            -m*(yg*q+zg*r) m*(yg*p+w) m*(zg*p-v) 0 -Iyz*q-Ixz*p+Iz*r Iyz*r+Ixy*p-Iy*q;
            m*(xg*q-w) -m*(zg*r+xg*p) m*(zg*q+u) Iyz*q+Ixz*p-Iz*r 0 -Ixz*r-Ixy*q+Ix*p;
            m*(xg*r+v) m*(yg*r-u) m*(xg*p+yg*q) -Iyz*r-Ixy*p+Iy*q Ixz*r+Ixy *q-Ix*p 0];
    end

%% Nested added-mass centripetal and Coriolis matrix generator function
    function Ca = genCa(A,nu)
        Asym = 0.5*(A+A');
        Ca = [zeros(3,3) -Smtrx(Asym(1:3,1:3)*nu(1:3,1)+Asym(1:3,4:6)*nu(4:6,1));
            -Smtrx(Asym(1:3,1:3)*nu(1:3,1)+Asym(1:3,4:6)*nu(4:6,1)) -Smtrx(Asym(4:6,1:3)*nu(1:3,1)+Asym(4:6,4:6)*nu(4:6,1))];
    end

%% Function body
Crb1 = genCrb(Mrb(1:6,1:6),nu(1:6,1),rg1);
Crb2 = genCrb(Mrb(7:12,7:12),nu(7:12,1),rg2);
Crb = [Crb1 zeros(6,6);zeros(6,6) Crb2];

Ca11 = genCa(Ainf(1:6,1:6),nu(1:6,1));
Ca12 = genCa(Ainf(1:6,7:12),nu(7:12,1));
Ca21 = genCa(Ainf(7:12,1:6),nu(1:6,1));
Ca22 = genCa(Ainf(7:12,7:12),nu(7:12,1));

Ca = [Ca11 Ca12;Ca21 Ca22];
M = Mrb + Ainf;
D = Ca;

Jth = zeros(6,6,2);
for k1 = 1:2
    phi = eta(4+(k1-1)*6,1);
    theta = eta(5+(k1-1)*6,1);
    psi = eta(6+(k1-1)*6,1);
    
    Rth = [cos(psi)*cos(theta) -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi) sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta);
        sin(psi)*cos(theta) cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi) -cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi);
        -sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi)];
    
    Tth = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
        0 cos(phi) -sin(phi);
        0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
    
    Jth(:,:,k1) = [Rth zeros(3);zeros(3) Tth];
end
eta1p = Jth(:,:,1)*nu(1:6,1);
eta2p = Jth(:,:,2)*nu(7:12,1);
etap = [eta1p;eta2p];


% Geta = [Jth(:,:,1)\G(1:6,1:6)*eta(1:6);Jth(:,:,2)\G(7:12,7:12)*eta(7:12)];
nup = M\(-(Crb+D)*nu-mu-G*eta+tau);
end