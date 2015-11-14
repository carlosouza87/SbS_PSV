function [Xs,Ys,Ns] = hydrint_brix93(Cxs,Cys,Cns,Lm,Tm,um,rho,ksi)

Xsmax = Cxs*(rho/2)*um^2*Lm*Tm;
Ysmax = Cys*(rho/2)*um^2*Lm*Tm;
Nsmax = Cns*(rho/2)*um^2*Lm*Tm;

if abs(ksi) <= 1
    ksi_values = -1.000:0.250:1.000;
    relXs = [-0.289 -0.690 -1.000 -0.850 -0.250 0.590 0.980 0.810 0.330];
    relYs = [0.289 0.345 -0.060 -0.595 -0.935 -0.982 -0.637 -0.250 -0.089];
    relNs = [0.264 0.706 1.000 0.873 0.221 -0.682 -0.927 -0.706 -0.424];

    Xs_Xsmax = interp1(ksi_values,relXs,ksi);
    Ys_Ysmax = interp1(ksi_values,relYs,ksi);
    Ns_Nsmax = interp1(ksi_values,relNs,ksi);

    Xs = Xsmax*Xs_Xsmax;
    Ys = Ysmax*Ys_Ysmax;
    Ns = Nsmax*Ns_Nsmax;
else
    Xs = 0;
    Ys = 0;
    Ns = 0;
end