function [dy,y,nsys,nn] = eqsim(y,t,newdt,nn)

global ktime data variable
global iwaves1st iwavesmd iwind icontrol_psv ihdsuction imemory icoupling icurr imooring

%% General constants
rho_water = data.const.rho_water;
rho_air = data.const.rho_air;
nu_water = data.const.nu_water;
g = data.const.g;

d2r = pi/180;       % Conversion factor, degrees to radians [rad/deg]
r2d = 180/pi;       % Conversion factor, radians to degrees [deg/rad]
knt2ms = 0.5144;    % Conversion factor, knots to meters per second [m/knt.s]
ms2knt = 1/knt2ms;  % Conversion factor, meters per second to knots [knt.s/m]

%% Ships properties
% Main dimensions
neq1=24;
neq2=18;

Loa(1) = data.ship(1).Loa;
Lpp(1) = data.ship(1).Lpp;
B(1) = data.ship(1).B;
T(1) = data.ship(1).T;
S(1) = data.ship(1).S;
Cb(1) = data.ship(1).Cb;
Cy(1) = data.ship(1).Cy;
lL(1) = data.ship(1).lL;
rg(:,1) = data.ship(1).rg;
Loa(2) = data.ship(2).Loa;
Lpp(2) = data.ship(2).Lpp;
B(2) = data.ship(2).B;
T(2) = data.ship(2).T;
S(2) = data.ship(2).S;
Cb(2) = data.ship(2).Cb;
Cy(2) = data.ship(2).Cy;
lL(2) = data.ship(2).lL;
rg(:,2) = data.ship(2).rg;

% Rigid body inertia matrices
Mrb1 = data.ship(1).Mrb;
Mrb2 = data.ship(2).Mrb;

% Restoration matrices
Ghd1 = data.ship(1).Ghd;
Ghd2 = data.ship(2).Ghd;


%% State vector
ypos = 0; % Current index of y vector
if idof == 1
    ndof = 6;
elseif idof == 2
    ndof = 3;
end

eta = y(1:2*ndof,1); % Positions
nu = y(2*ndof+1:4*ndof,1); % Velocities
eta1 = eta(1:ndof); % Vessel 1 (FPSO) positions
eta2 = eta(ndof+1:2*ndof); % Vessel 2 (PSV) positions
nu1 = nu(1:ndof);  % Vessel 1 (FPSO) velocities
nu2 = nu(ndof+1:2*ndof); % Vessel 2 (PSV) velocities

ypos = ypos + 4*ndof; % Update ypos

%% FPSO mooring system

if imooring == 1
    load moor_coord.txt % File with anchors and fairlead coordinates
    load moor_lxT.txt % File with horizontal distance vs. line tensions
    
    % Vector with current CG horizontal position
    if idof == 1
        eta_hor = [eta1(1);eta1(2);eta1(6)];
    elseif idof == 2
        eta_hor = [eta1(1);eta1(2);eta1(3)];
    end
    [Xmoor,Ymoor,Nmoor,broke] = mooring(moor_coord,moor_lxT,eta_hor);
        
    tau_moor = [Xmoor;Ymoor;Nmoor]; % Vector with mooring loads
    
    if sum(broke) > 0
        lines_brk = 1:length(broke);
        lines_brk = lines_brk(broke==1);
        warning_msg = ['Broken lines (' num2str(lines_brk) ')'];
        warning(warning_msg);
    end
        
    
    variable.tau_amarras(:,ktime) = tau_amarras;
    
else
    tau_amarras = zeros(12,1);
    
end

%% Particular loads
tau_waves1st = zeros(12,1);
tau_wavesmd = zeros(12,1);
tau_wind = zeros(12,1);
tau_curr = zeros(12,1);
tau_ctr = zeros(12,1);
tau_hdsuction = zeros(12,1);
for k1 = 1:2
    
    u = nu(1+(k1-1)*6,1);
    v = nu(2+(k1-1)*6,1);
    psi = eta(6+(k1-1)*6,1);
    
    % Waves - 1st order
    if iwaves1st == 1
        betaw = d2r*data.environment.betaw;
        betaFTF = data.ship(k1).waveincid;
        Fwf = data.ship(k1).Fwf;
        beta = 180/pi*norm02pi(betaw-psi); % wave incidence direction [deg]
        for k2 = 1:6
            if (k2 == 1) && (t >= 2.6)
                debug = 1;
            end
            tau_waves1st(k2+6*(k1-1),1) = interp1(betaFTF,Fwf(ktime,:,k2),beta,'linear','extrap');
        end
    end
    variable.ship(k1).waves1st(:,ktime) = tau_waves1st(1+(k1-1)*6:6+(k1-1)*6);
    
    %Waves - mean drift
    if iwavesmd == 1
        u_md = variable.ship(k1).u_md;
        v_md = variable.ship(k1).v_md;
        psi_md = variable.ship(k1).psi_md;
        if abs(psi_md-psi) > 0.0873
            variable.ship(k1).calc_mdrift = 1;
            variable.ship(k1).psi_md = psi;
        elseif sqrt((u_md-u)^2+(v_md-v)^2) > 0.3
            variable.ship(k1).calc_mdrift = 1;
            variable.ship(k1).u_md = nu(1);
            variable.ship(k1).v_md = nu(2);
        else
            variable.ship(k1).calc_mdrift = 0;
        end
        if variable.ship(k1).calc_mdrift == 1
            betaw = d2r*data.environment.betaw;
            betaWMD = data.ship(k1).waveincid; %1)olhar tamanho 2) fazero caminho17 vetor
            Fmd = data.ship(k1).Fmd;
            beta = 180/pi*norm02pi(betaw-psi); % wave incidence direction [deg]
            if beta <= 180
                tau_wavesmd(1+6*(k1-1),1) = interp1(betaWMD,Fmd(:,1),beta,'linear');
                tau_wavesmd(2+6*(k1-1),1) = interp1(betaWMD,Fmd(:,2),beta,'linear');
                tau_wavesmd(6+6*(k1-1),1) = interp1(betaWMD,Fmd(:,3),beta,'linear');
            else
                tau_wavesmd(1+6*(k1-1),1) = interp1(betaWMD,Fmd(:,1),beta,'linear'); %arquivo data.ship(k1).waveincid;360 -> data.ship(k1).waveincid;360 - betaWMD
                tau_wavesmd(2+6*(k1-1),1) = -interp1(betaWMD,Fmd(:,2),beta,'linear');
                tau_wavesmd(6+6*(k1-1),1) = -interp1(betaWMD,Fmd(:,3),beta,'linear');
            end
            variable.ship(k1).lastktime_md = ktime;
        else
            lastktime_md = variable.ship(k1).lastktime_md;
            tau_wavesmd(1+6*(k1-1):6*k1,1) = variable.ship(k1).tau_wavesmd(:,lastktime_md);
        end
    else
        tau_wavesmd = zeros(12,1);
    end
    variable.ship(k1).tau_wavesmd(:,ktime) = tau_wavesmd(1+6*(k1-1):6*k1,1);
    
    
    % Wind
    if iwind == 1
        u_wn = variable.ship(k1).u_wn;
        v_wn = variable.ship(k1).v_wn;
        psi_wn = variable.ship(k1).psi_wn;
        if abs(psi_wn-psi) > 0.0873
            variable.ship(k1).calc_wind = 1;
            variable.ship(k1).psi_wn = psi;
        elseif sqrt((u_wn-u)^2+(v_wn-v)^2) > 0.3
            variable.ship(k1).calc_wind = 1;
            variable.ship(k1).u_wn = nu(1);
            variable.ship(k1).v_wn = nu(2);
        else
            variable.ship(k1).calc_wind = 0;
        end
        if variable.ship(k1).calc_wind == 1
            bow = data.ship(k1).bow;
            load = data.ship(k1).load;
            At = data.ship(k1).At;
            Al = data.ship(k1).Al;
            coefwind = data.environment.coefwind;
            gammaw = d2r*data.environment.gammaw;
            Uwind = data.environment.Uw;
            gamma = gammaw - psi;
            
            uw = Uwind*cos(gamma) - u;
            vw = Uwind*sin(gamma) - v;
            
            gamma = atan2(-vw,-uw);
            
            % Angle normalization
            if abs(gamma) > 2*pi
                gamma = rem(gamma,2*pi);
            end
            if gamma < 0
                gamma = 2*pi + gamma;
            end
            Uw_r = sqrt(uw^2+vw^2); % wind velocity related to the ship hull
            [Xwind,Ywind,Nwind] = ocimf77(coefwind,r2d*gamma,rho_air,Uw_r,At,Al,Lpp(k1),bow,load);    % wind forces and moment [N, N, Nm]
            tau_wind(1+6*(k1-1),1) = Xwind;
            tau_wind(2+6*(k1-1),1) = Ywind;
            tau_wind(6+6*(k1-1),1) = Nwind;
            variable.ship(k1).lastktime_wn = ktime;
        else
            lastktime_wn = variable.ship(k1).lastktime_wn;
            tau_wind = variable.ship(k1).tau_wind(:,lastktime_wn);
        end
    end
    variable.ship(k1).tau_wind(:,ktime) = tau_wind;
    
    
    
    % Current
    if icurr == 1
        %     if k1 == 1
        alphac = d2r*data.environment.alphac;
        alpha = norm02pi(alphac-psi);   % incidence direction, normalized between 0 and 2*pi [rad]
        Uc = data.environment.Uc;
        ur = u - Uc*cos(alpha);
        vr = v - Uc*sin(alpha);
        [Xcurr,Ycurr,Ncurr] = current(Lpp(k1),B(k1),T(k1),S(k1),Cb(k1),Cy(k1),lL(k1),ur,vr,nu(6+(k1-1)*6));
        tau_curr(1+6*(k1-1),1) = Xcurr;
        tau_curr(2+6*(k1-1),1) = Ycurr;
        tau_curr(6+6*(k1-1),1) = Ncurr;
        %     else
        %         alphac = d2r*data.environment.alphac;
        %         alpha = norm02pi(alphac-psi);   % incidence direction, normalized between 0 and 2*pi [rad]
        %         Uc = data.environment.Uc;
        %
        %
        %         [Xcurr,Ycurr,Ncurr] = correnteza_ipt(Lpp(k1),T(k1),Uc,alphac,rho_water);
        %         tau_curr(1+6*(k1-1),1) = Xcurr;
        %         tau_curr(2+6*(k1-1),1) = Ycurr;
        %         tau_curr(6+6*(k1-1),1) = Ncurr;
        %
        %     end
        
    end
    
    variable.ship(k1).tau_curr(:,ktime) = tau_curr;
    
    
    %% Control
    if icontrol_psv == 1
        
        y_m=[eta2(1); eta2(2); eta2(6)];
        variable.eta_hat(:,ktime) = [0;0;0;0;0;0;eta_hat(1,1);eta_hat(2,1);0;0;0;eta_hat(3,1)];
        variable.nu_hat(:,ktime) = [0;0;0;0;0;0;nu_hat(1,1);nu_hat(2,1);0;0;0;nu_hat(3,1)];
        
        % PSV control parameters
        eta_hat(1,1) = y(ypos+1,1);
        eta_hat(2,1) = y(ypos+2,1);
        eta_hat(3,1) = y(ypos+3,1);
        nu_hat(1,1) = y(ypos+4,1);
        nu_hat(2,1) = y(ypos+5,1);
        nu_hat(3,1) = y(ypos+6,1);
        ksi_hat(1,1) = y(ypos+7,1);
        ksi_hat(2,1) = y(ypos+8,1);
        ksi_hat(3,1) = y(ypos+9,1);
        ksi_hat(4,1) = y(ypos+10,1);
        ksi_hat(5,1) = y(ypos+11,1);
        ksi_hat(6,1) = y(ypos+12,1);
        b_hat(1,1) = y(ypos+13,1);
        b_hat(2,1) = y(ypos+14,1);
        b_hat(3,1) = y(ypos+15,1);
        variable.ship(2).eta(:,1) = y(7:12,1);
        
        eta_ref=[data.ship(2).control.xref; data.ship(2).control.yref;data.ship(2).control.psiref];
        nu_ref=[0;0;0];
        
        eta_error= eta_ref-eta_hat;
        nu_error=nu_ref-nu_hat;
        integra_eta_error=[y(neq1+neq2-2,1);y(neq1+neq2-1,1);y(neq1+neq2,1)];
        
        
        psi=eta2(6,1);
        Jpsi = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];% rotation matrix helio
        
        
        Kp_psv_x = data.ship(2).control.Kp_x;
        Kd_psv_x = data.ship(2).control.Kd_x;
        Ki_psv_x = data.ship(2).control.Ki_x;
        
        Kp_psv_y = data.ship(2).control.Kp_y;
        Kd_psv_y = data.ship(2).control.Kd_y;
        Ki_psv_y = data.ship(2).control.Ki_y;
        
        Kp_psv_psi = data.ship(2).control.Kp_psi;
        Kd_psv_psi = data.ship(2).control.Kd_psi;
        Ki_psv_psi = data.ship(2).control.Ki_psi;
        
        Kp_psv=diag([Kp_psv_x;Kp_psv_y;Kp_psv_psi]);
        Kd_psv=diag([Kd_psv_x;Kd_psv_y;Kd_psv_psi]);
        Ki_psv=diag([Ki_psv_x;Ki_psv_y;Ki_psv_psi]);
        
        tau_Kp=Kp_psv*eta_error;
        tau_Ki=Ki_psv*integra_eta_error;
        tau_Kd=Kd_psv*nu_error;
        
        
        b_hat_movel=inv(Jpsi)*b_hat;
        tau_cont=inv(Jpsi)*tau_Kp+inv(Jpsi)*tau_Ki+tau_Kd-b_hat_movel;
        tau_cont_aux=inv(Jpsi)*tau_Kp+inv(Jpsi)*tau_Ki+tau_Kd;
        tau_cont_aux_Kp=inv(Jpsi)*tau_Kp;
        tau_cont_aux_Kd=tau_Kd;
        
        % Control vector
        tau_ctr(7:12,1) = [tau_cont(1,1);tau_cont(2,1);0;0;0;tau_cont(3,1)];
        
    else
        eta_error= [0;0;0];
        nu_error=[0;0;0];
        integra_eta_error=[0;0;0];
        Jpsi = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
        b_hat_movel=inv(Jpsi)*b_hat;
        tau_cont=[0;0;0];
        tau_ctr(7:12,1) = [tau_cont(1,1);tau_cont(2,1);0;0;0;tau_cont(3,1)];
        tau_cont_aux=[0;0;0];
        tau_cont_aux_Kp=[0;0;0];
        tau_cont_aux_Kd=[0;0;0];
    end
    
    
    %% Suction interaction loads
    
    
    %CARLOS - FORMULACAO UTILIZANDO VANTORRE
    
    if ihdsuction == 1
        
        alphac = d2r*data.environment.alphac;
        alpha = norm02pi(alphac-psi);   % incidence direction, normalized between 0 and 2*pi [rad]
        Uc = data.environment.Uc;
        
        if k1 == 1
            xcc = (eta1(1)-eta2(1))*cos(eta1(6)) + (eta1(2)-eta2(2))*sin(eta1(6));
            ksi = xcc/(0.5*(Lpp(1)+Lpp(2)));
            
            U1 = nu(1) - Uc*cos(norm02pi(alphac - eta(6,1)));
            U2 = nu(7) - Uc*cos(norm02pi(alphac - eta(12,1)));
            
            [Xh,Yh,Nh] = hydrint_vantorre01(ksi,rho_water,Lpp(1),B(1),T(1),U1,U2);
            tau_hdsuction(1:6,1) = [0;Yh;0;0;0;Nh];
        elseif k1 == 2
            [Xh,Yh,Nh] = hydrint_vantorre01(ksi,rho_water,Lpp(2),B(2),T(2),U1,U2);
            Yh = -Yh;
            Nh = -Nh;
            tau_hdsuction(7:12,1) = [0;Yh;0;0;0;Nh];
        end
    end
    variable.ship(1).tau_hdsuction(:,ktime) = tau_hdsuction(1:6,1);
    variable.ship(2).tau_hdsuction(:,ktime) = tau_hdsuction(7:12,1); %tau_hdsuction
    
end
nn = nn + 1;

%% Observer parameters
zetanotch1 = 1; zetanotch2 = 1; zetanotch3 = 1;
zeta1 = 0.1; zeta2 = 0.1; zeta3 = 0.1;

Tpf1=9.0;%9.5;
Tpf2=9.0;
Tpf3=9.0; %10%10.5;%9.5%10.5;
Tb=diag([500,500,500]);

w01 = 2*pi/(Tpf1); w02 = 2*pi/(Tpf2); w03 = 2*pi/(Tpf3);
wc1 = 1.1*w01; wc2 = 1.1*w02; wc3 = 1.1*w03;%1.1*w03;

%equacoes (4.23 a 4.26)
Ka1(1,1) = -2*(zetanotch1-zeta1)*wc1/w01;
Ka1(2,2) = -2*(zetanotch2-zeta2)*wc2/w02;
Ka1(3,3) = -2*(zetanotch3-zeta3)*wc3/w03;
Ka1(4,1) = 2*(zetanotch1-zeta1)*w01;
Ka1(5,2) = 2*(zetanotch2-zeta2)*w02;
Ka1(6,3) = 2*(zetanotch3-zeta3)*w03;
Ka2(1,1) = wc1;
Ka2(2,2) = wc2;
Ka2(3,3) = wc3;
Ka3=1e5*[1 0 0; 0 1 0; 0 0 1e4];%*1e3;
Ka4=Ka3*100;

OMEGA21 = -diag([w01^2;w02^2;w03^2]); OMEGA22 = -diag([2*zeta1*w01;2*zeta2*w02;2*zeta3*w03]);
OMEGA = [zeros(3) eye(3);OMEGA21 OMEGA22];
GAMMA = [zeros(3) eye(3)];



%% Equations of motions

if isimtype == 1
    % Convolution integral calculation for Cummins equation
    
    % Evaluation of retardation functions and infinite added mass
    
    K11 = variable.K11;
    K12 = variable.K12;
    K21 = variable.K21;
    K22 = variable.K22;
    
    Tij11=data.hydro.Tij11;
    Tij12=data.hydro.Tij12;
    Tij21=data.hydro.Tij21;
    Tij22=data.hydro.Tij22;
    
    size_k = [size(K11,3) size(K12,3) size(K21,3) size(K22,3)];
    size_k = max (size_k);
    
    K11 = zeros(6,6,size_k);
    K12 = zeros(6,6,size_k);
    K21 = zeros(6,6,size_k);
    K22 = zeros(6,6,size_k);
    
    
    K11(:,:,1:size(variable.K11,3)) = variable.K11(:,:,1:size(variable.K11,3));
    K12(:,:,1:size(variable.K12,3)) = variable.K12(:,:,1:size(variable.K12,3));
    K21(:,:,1:size(variable.K21,3)) = variable.K21(:,:,1:size(variable.K21,3));
    K22(:,:,1:size(variable.K22,3)) = variable.K22(:,:,1:size(variable.K22,3));
    
    
    K=[K11 K12 ; K21 K22];
    
    Tij_max = [Tij11 Tij12 ; Tij21 Tij22];
    Tij_max = Tij_max';
    
    Tij_max = max(Tij_max);
    
    [mu] = convolution_integral(K,Tij_max,ktime,nu,t,newdt);
end

tau_b_hat=[0;0;0;0;0;0;b_hat_movel(1,1);b_hat_movel(2,1);0;0; 0; b_hat_movel(3,1)];

tau_ext = tau_waves1st + tau_wavesmd + tau_wind + tau_curr + tau_hdsuction + tau_amarras;
tau = tau_ext + tau_ctr;

variable.tau_ctr(:,ktime) = tau_ctr;
variable.tau_b_hat_psv(:,ktime)=tau_b_hat;
variable.tau_ext(:,ktime)=tau_ext;
variable.tau_cont_aux(:,ktime)=tau_cont_aux;
variable.tau_cont_aux_Kp(:,ktime)=tau_cont_aux_Kp;
variable.tau_cont_aux_Kd(:,ktime)=tau_cont_aux_Kd;


psi=eta2(6);
Jpsi = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];% rotation matrix helio

tau_ctr_alt=[tau_ctr(7);tau_ctr(8);tau_ctr(12)]; %matriz criada que deve ser 3x1
MRB2_alt=[MRB2(1,1)*1.1,MRB2(1,2),MRB2(1,6);MRB2(2,1),MRB2(2,2)*1.9,MRB2(2,6);MRB2(6,1),MRB2(6,2),MRB2(6,6)*1.2];%colocando massas adicionais em x, y e psi

[etap_hat,nup_hat,ksip_hat,bp_hat,eps] = nonl_observer(y_m,Jpsi,eta_hat,nu_hat,ksi_hat,b_hat,MRB2_alt,D,tau_ctr_alt,OMEGA,GAMMA,Ka1,Ka2,Ka3,Ka4,Tb,t);

MRB = [MRB1 zeros(6,6);zeros(6,6) MRB2];
G = [G1 zeros(6,6);zeros(6,6) G2];

% Aqui, trocar imemory por isimtype, eqmot_memory por eqmot_cummins e calcular o mu (Carlos, 08/12/15)
if imemory ==1
    Ainf = data.hydro.A_inf;
    [etap,nup] = eqmot_memory(eta,nu,rg(:,1),rg(:,2),MRB,Ainf,G,mu,tau,newdt);
else
    A_fixfreq = [data.hydro.A11_fixfreq data.hydro.A12_fixfreq;data.hydro.A21_fixfreq data.hydro.A22_fixfreq]; %espa�o nova coluna e ; nova linha
    B_fixfreq = [data.hydro.B11_fixfreq data.hydro.B12_fixfreq;data.hydro.B21_fixfreq data.hydro.B22_fixfreq];
    nu_m = [0*5*0.5144;0;0;0;0;0;0*5*0.5144;0;0;0;0;0]; %foram zerados os primeiro e o s�timo termo pois eles s�o a vel. avan�o que queremos nula.
    [etap,nup] = eqmot_fixfreq(eta,nu,nu_m,rg(:,1),rg(:,2),MRB,A_fixfreq,B_fixfreq,G,tau); %linha que foi zerada (manter zerada muito problem�tica)
end

variable.etap(:,ktime) = etap(:,1);
variable.nu(:,ktime) = nu(:,1);
dt=variable.dt;

% if newdt==1
% i=5;
% figure(35)
% t_plot= (1:size(variable.nu,2))*dt;
% t_plot=t_plot';
% % K_plot(:,1)= K(i,i,:);
% plot(t_plot,variable.nu(i,1:length(t_plot)),'o-','LineWidth',2)
% xlabel ('tempo')
% ylabel (['variable.nu(' num2str(i) ',' num2str(i) ')'])
% grid on
% end

nsys = size(y,1);
dy = zeros(nsys,1);
ypos = 0;

dy(1:12,1) = etap;
dy(13:24,1) = nup;

ypos=24;
for iandrey = 1:3
    dy(ypos+iandrey,1)=etap_hat(iandrey);
    dy(ypos+3+iandrey,1)=nup_hat(iandrey);
    dy(ypos+6+iandrey,1)=ksip_hat(iandrey);
    dy(ypos+9+iandrey,1)=ksip_hat(iandrey+3);
    dy(ypos+12+iandrey,1)=bp_hat(iandrey);
    dy(ypos+15+iandrey,1)=eta_error(iandrey);
end


ypos = neq1+neq2;  %ypos + 24+15 Helio

