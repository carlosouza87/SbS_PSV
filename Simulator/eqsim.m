function [dy,y,nsys] = eqsim(y,t,newdt)

global ktime data variable flag

% Read flags
isimtype = flag.isimtype;
idof =  flag.idof;
iwaves1st = flag.iwaves1st;
iwavesmd = flag.iwavesmd;
iwind = flag.iwind;
icurr = flag.icurr;
icontrol_psv = flag.icontrol_psv;
ihdsuction = flag.ihdsuction;
imooring = flag.imooring;

%% General constants
rho_water = data.constants.rho_water;
rho_air = data.constants.rho_air;
nu_water = data.constants.nu_water;
g = data.constants.g;

d2r = pi/180;       % Conversion factor, degrees to radians [rad/deg]
r2d = 180/pi;       % Conversion factor, radians to degrees [deg/rad]
knt2ms = 0.5144;    % Conversion factor, knots to meters per second [m/knt.s]
ms2knt = 1/knt2ms;  % Conversion factor, meters per second to knots [knt.s/m]

% Number of degrees of freedom
if idof == 1
    ndof = 6;
elseif idof == 2
    ndof = 3;
end

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
    
    
    variable.tau_moor(:,ktime) = tau_moor;
    
else
    tau_moor = zeros(ndof*2,1);
end

%% Particular loads
tau_waves1st = zeros(ndof*2,1);
tau_wavesmd = zeros(ndof*2,1);
tau_wind = zeros(ndof*2,1);
tau_curr = zeros(ndof*2,1);
tau_ctr = zeros(ndof*2,1);
tau_hdsuction = zeros(ndof*2,1);
for k1 = 1:2
    
    u = nu(1+(k1-1)*ndof,1);
    v = nu(2+(k1-1)*ndof,1);
    psi = eta(ndof+(k1-1)*ndof,1);
    
    % Waves - 1st order
    if iwaves1st == 1
        betaw = d2r*data.environment.betaw;
        betaw1st = data.ship(k1).waves_incid;
        Fwf = data.ship(k1).Fwf;
        beta = 180/pi*norm02pi(betaw-psi); % wave incidence direction [deg]
        for k2 = 1:ndof
            tau_waves1st(k2+ndof*(k1-1),1) = interp1(betaw1st,Fwf(ktime,:,k2),beta,'linear','extrap');
        end
    end
    variable.ship(k1).waves1st(:,ktime) = tau_waves1st(1+(k1-1)*ndof:ndof+(k1-1)*ndof);
    
    % Waves - mean drift
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
            betaWMD = data.ship(k1).waves_incid;
            Fmd = data.ship(k1).Fmd;
            beta = 180/pi*norm02pi(betaw-psi); % wave incidence direction [deg]
            if beta <= 180
                tau_wavesmd(1+ndof*(k1-1),1) = interp1(betaWMD,Fmd(:,1),beta,'linear');
                tau_wavesmd(2+ndof*(k1-1),1) = interp1(betaWMD,Fmd(:,2),beta,'linear');
                tau_wavesmd(ndof+ndof*(k1-1),1) = interp1(betaWMD,Fmd(:,3),beta,'linear');
            else
                tau_wavesmd(1+ndof*(k1-1),1) = interp1(betaWMD,Fmd(:,1),beta,'linear'); %arquivo data.ship(k1).waveincid;360 -> data.ship(k1).waveincid;360 - betaWMD
                tau_wavesmd(2+ndof*(k1-1),1) = -interp1(betaWMD,Fmd(:,2),beta,'linear');
                tau_wavesmd(ndof+ndof*(k1-1),1) = -interp1(betaWMD,Fmd(:,3),beta,'linear');
            end
            variable.ship(k1).lastktime_md = ktime;
        else
            lastktime_md = variable.ship(k1).lastktime_md;
            tau_wavesmd(1+ndof*(k1-1):ndof*k1,1) = variable.ship(k1).tau_wavesmd(:,lastktime_md);
        end
    else
        tau_wavesmd = zeros(2*ndof,1);
    end
    variable.ship(k1).tau_wavesmd(:,ktime) = tau_wavesmd(1+ndof*(k1-1):ndof*k1,1);
    
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
            load_cond = data.ship(k1).load;
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
            [Xwind,Ywind,Nwind] = ocimf77(coefwind,r2d*gamma,rho_air,Uw_r,At,Al,Lpp(k1),bow,load_cond);    % wind forces and moment [N, N, Nm]
            tau_wind(1+ndof*(k1-1),1) = Xwind;
            tau_wind(2+ndof*(k1-1),1) = Ywind;
            tau_wind(ndof+ndof*(k1-1),1) = Nwind;
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
        [Xcurr,Ycurr,Ncurr] = current(Lpp(k1),B(k1),T(k1),S(k1),Cb(k1),Cy(k1),lL(k1),ur,vr,nu(ndof+(k1-1)*ndof));
        tau_curr(1+ndof*(k1-1),1) = Xcurr;
        tau_curr(2+ndof*(k1-1),1) = Ycurr;
        tau_curr(ndof+ndof*(k1-1),1) = Ncurr;
        %     else
        %         alphac = d2r*data.environment.alphac;
        %         alpha = norm02pi(alphac-psi);   % incidence direction, normalized between 0 and 2*pi [rad]
        %         Uc = data.environment.Uc;
        %
        %
        %         [Xcurr,Ycurr,Ncurr] = correnteza_ipt(Lpp(k1),T(k1),Uc,alphac,rho_water);
        %         tau_curr(1+ndof*(k1-1),1) = Xcurr;
        %         tau_curr(2+ndof*(k1-1),1) = Ycurr;
        %         tau_curr(ndof+ndof*(k1-1),1) = Ncurr;
        %
        %     end
        
    end
    
    variable.ship(k1).tau_curr(:,ktime) = tau_curr;
    
    
    %% Control
    if icontrol_psv == 1
        
        y_m=[eta2(1); eta2(2); eta2(ndof)];
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
        variable.ship(2).eta(:,1) = y(ndof+1:2*ndof,1);
        
        % PSV observer parameters
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
        
        eta_ref = [data.ship(2).control.xref; data.ship(2).control.yref;data.ship(2).control.psiref];
        nu_ref = [0;0;0];
        
        eta_error= eta_ref-eta_hat;
        nu_error=nu_ref-nu_hat;
        integra_eta_error=[y(neq1+neq2-2,1);y(neq1+neq2-1,1);y(neq1+neq2,1)];
        
        
        psi = eta2(ndof,1);
        Jpsi = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
        
        
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
        
        tau_b_hat=[0;0;0;0;0;0;b_hat_movel(1,1);b_hat_movel(2,1);0;0; 0; b_hat_movel(3,1)];
        
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
        
    else
        %         eta_error= [0;0;0];
        %         nu_error=[0;0;0];
        %         integra_eta_error=[0;0;0];
        %         Jpsi = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
        %         b_hat_movel=inv(Jpsi)*b_hat;
        %         tau_cont=[0;0;0];
        %         tau_ctr(ndof+1:2*ndof,1) = [tau_cont(1,1);tau_cont(2,1);zeros(ndof-3,1);tau_cont(3,1)];
        %         tau_cont_aux=[0;0;0];
        %         tau_cont_aux_Kp=[0;0;0];
        %         tau_cont_aux_Kd=[0;0;0];
    end
    
    
    
end


%% Equations of motions
    % Rigid-body inertia matrix for both ships
    Mrb = [Mrb1 zeros(6,6);zeros(6,6) Mrb2];
    
    % Hydrostatic restoration matrix for both ships
    Ghd = [Ghd1 zeros(6,6);zeros(6,6) Ghd2];
if isimtype == 1
     
    % Retardation functions from previously generated structure
    K=[data.hydro.K11 data.hydro.K12; data.hydro.K21 data.hydro.K22];

    % Calculation of convolution integral for Cummins equation
    if newdt == 1
        nu_mem = [variable.ship(1).nu(:,1:ktime);variable.ship(2).nu(:,1:ktime)];
        [mu] = convolution_integral(K,nu_mem,t,ktime);
%      mu = zeros(12,1);
        variable.mu(:,ktime) = mu;
    else
        [mu] = variable.mu(:,ktime);
    end
    
    % Infinite-frequency added inertia matrix for both ships
    A_inf = [data.hydro.A11_inf data.hydro.A12_inf;
         data.hydro.A21_inf data.hydro.A22_inf];   
     
    % Vector of external loads
    tau_ext = tau_waves1st + tau_wavesmd + tau_wind + tau_curr + tau_hdsuction + tau_moor;
    tau = tau_ext + tau_ctr;  
    
    % Equations of motions
    [etap,nup] = eqmotions_6(eta,nu,rg(:,1),rg(:,2),Mrb,A_inf,Ghd,mu,tau);

elseif isimtype == 2
    % Vector mu is dummy for isimtype and all elements must be set to zero
    mu_dummy = zeros(2*ndof,1);
    
    % Zero-frequency added inertia matrix
    A_0 = [data.hydro.A11 data.hydro.A12;
         data.hydro.A21 data.hydro.A22]; 
     
    % Vector of external loads
    tau_ext = tau_wavesmd + tau_wind + tau_curr + tau_hdsuction + tau_moor;
    tau = tau_ext + tau_ctr;  
    
    % Equations of motions
    [etap,nup] = eqmotions_6(eta,nu,rg(:,1),rg(:,2),Mrb,A_0,Ghd,mu_dummy,tau);
end

variable.etap(:,ktime) = etap;
variable.nup(:,ktime) = nup;

nsys = size(y,1);
dy = zeros(nsys,1);
dypos = 0;
dy(1:2*ndof,1) = etap;
dypos = dypos + 2*ndof;
dy(dypos+1:2*dypos,1) = nup;



