function [Xmoor,Ymoor,Nmoor,broke] = mooring(coord,lxT,eta_hor)
% [Xmoor,Ymoor,Nmoor,broke] = mooring(coord,lxT,eta_hor)
% Calculate horizontal mooring loads based on fairleads and anchors
% coordinates, and on a lookup table with line lengths and the respective
% induced tensions. It is admitted that all lines have the same properties.
% Input is:
% coord - N x 4 matrix with anchor Earth-fixed coordinates (xa, ya)  
% and fairlead body-fixed coordinates (xf_bf, yf_bf). The coordinates
% should be disposed as follows:
% xa1 ya1 xf_bf1 yf_bf1
% xa2 ya2 xf_bf2 yf_bf2
%  .   .    .      .
%  .   .    .      .
%  .   .    .      .
% xaN yaN xf_bfN yf_bfN
% 
% lxT - M x 2 matrix with a lookup table relating the horizontal distance 
% from a fairlead to the respective anchor and the tension induced in the
% line. M is the number of distance x tension values available in the
% table. The distances and tensions should be disposed as follows:
% l1 T1
% l2 T2
% .  .
% .  .
% .  .
% lM TM
% 
% eta_hor - 3 x 1 vector with the vessel CG's horizontal posisitions (x and
% y) and heading (psi).
% 
% Output is:
% Xmoor, Ymoor and Nmoor - Mooring forces in X and Y and and moment around 
% Z axis in body-fixed coordinates.
% 
% broke - N x 1 vector with binary values indicating whether the
% corresponding line has broken (1) or not (0).


% CG coordinates and vessel heading
xg = eta_hor(1);
yg = eta_hor(2);
psi = eta_hor(3);

% Anchors coordinates
xa = coord(:,1);
ya = coord(:,2);

% Fairleads coordinates (in body-fixed frame)
xf_bf = coord(:,3);
yf_bf = coord(:,4);

% Number of lines
Nlines = size(coord,1);

% Lengths and tensions from lookup table
l = lxT(:,1);
T = lxT(:,2);

Xk_bf = zeros(1,Nlines);
Yk_bf = zeros(1,Nlines);
Nk_bf = zeros(1,Nlines);
broke = zeros(1,Nlines);

for k1=1:Nlines
    % Transform fairlead coordinates to Earth-fixed frame
    xf = xg + xf_bf(k1)*cos(psi)-yf_bf(k1)*sin(psi);
    yf = yg + xf_bf(k1)*sin(psi)+yf_bf(k1)*cos(psi);
    
    % Distance from current fairlead to respective anchor
    dist_fa = sqrt((xf-xa(k1))^2+(yf-ya(k1))^2);
    
    if dist_fa < min(l)
        T_line = 0;
    elseif dist_fa >= max(l)
        broke(k1) = 1;
        T_line = max(T);
    else
        T_line = interp1(l,T,dist_fa);
    end
    
    % Project line tension over X and Y directions (Earth-fixed frame)
    Xk = T_line*(xa(k1)-xf)/dist_fa; 
    Yk = T_line*(ya(k1)-yf)/dist_fa;
    
    % Transform loads to body-fixed frame
    Xk_bf(k1) = Xk*cos(psi) + Yk*sin(psi);
    Yk_bf(k1) = Yk*cos(psi) - Xk*sin(psi);
    
    % Moment
    Nk_bf(k1) = -Xk_bf(k1)*yf_bf(k1) + Yk_bf(k1)*xf_bf(k1);
end

Xmoor = sum(Xk_bf);
Ymoor = sum(Yk_bf);
Nmoor = sum(Nk_bf);


