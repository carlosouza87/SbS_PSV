function [etap,nup] = eqmot_cummins(eta,nu,rg1,rg2,Mrb,Ainf,G,mu,tau)
% function [etap,nup] = eqmot_cummins(eta,nu,rg1,rg2,Mrb,Ainf,G,mu,tau)
% Performs calculation of the equations of motion for both ships, with
% consideration of fluid memory effects. Inputs are:
%
% eta - Vector of positions for both vessels, earth-fixed reference system [12x1]
% nu - Vector of velocities for both vessels, ship-fixed reference system [12x1]
% rg1 - Position of center of gravity (ship 1) [1x3]
% rg2 - Position of center of gravity (ship 2) [1x3]
% Mrb - Rigid body inertia matrix for both vessels [12x12]
% Ainf - Infinity frequency added inertia matrix for both vessels, including couplings [12x12]
% G - Hydrostatic restoration matrix for both vessels [12x12]
% mu - Convolution integral of Cummins equation, calculated in an external function [12x1]
% tau - Vector with external forces and moment acting on both ships [12x1]
%
% Outputs are:
%
% etap - vector with time derivatives of eta
% nup - vector with time derivatives of nu
%
% Carlos Souza, 08/12/2015 - Universidade de São Paulo

%% Nested rigid-body centripetal and Coriolis matrix generator function
% See eq. 3.66 in (Fossen, 2002)
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
% See eq. 3.100 in (Fossen, 2002)
    function Ca = genCa(A,nu)
        Asym = 0.5*(A+A');
        Ca = [zeros(3,3) -Smtrx(Asym(1:3,1:3)*nu(1:3,1)+Asym(1:3,4:6)*nu(4:6,1));
            -Smtrx(Asym(1:3,1:3)*nu(1:3,1)+Asym(1:3,4:6)*nu(4:6,1)) -Smtrx(Asym(4:6,1:3)*nu(1:3,1)+Asym(4:6,4:6)*nu(4:6,1))];
    end

%% Function body
% Calculate rigid-body centripetal and Coriolis matrix for each vessel,
% and include both in a 12 x 12 matrix.
Crb1 = genCrb(Mrb(1:6,1:6),nu(1:6,1),rg1);
Crb2 = genCrb(Mrb(7:12,7:12),nu(7:12,1),rg2);
Crb = [Crb1 zeros(6,6);zeros(6,6) Crb2];

% Calculate added-mass centripetal and Coriolis matrix, using the
% submatrices of Ainf
Ca11 = genCa(Ainf(1:6,1:6),nu(1:6,1));
Ca12 = genCa(Ainf(1:6,7:12),nu(7:12,1));
Ca21 = genCa(Ainf(7:12,1:6),nu(1:6,1));
Ca22 = genCa(Ainf(7:12,7:12),nu(7:12,1));
Ca = [Ca11 Ca12;Ca21 Ca22];

M = Mrb + Ainf; % Total inertia matrix

% Matrix for transformation of coordinates from body-fixed to Earth-fixed frames
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

% Calculate etap from the current value of nu 
eta1p = Jth(:,:,1)*nu(1:6,1);
eta2p = Jth(:,:,2)*nu(7:12,1);
etap = [eta1p;eta2p];

% Equations of motions for calculating nup
nup = M\(-(Crb+Ca)*nu-mu-G*eta+tau);
end

% Reference:
% Fossen, T. I. (2002) - "Marine control systems - Guidance, Navigation and
% Control of Ships, Rigs and Underwater Vehicles". Marine Cybernetics,
% Trondheim, Norway.