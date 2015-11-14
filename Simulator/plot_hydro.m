clear all;close all;clc

load hydrodata_cs1d3m_D_3E5

w = data.hydro.freqs;
lw = length(w);

% dof = [1,1;1,2;1,3;1,4;1,5;1,6;2,1;2,2;2,3;2,4;2,5;2,6;3,1;3,2;3,3;3,4;3,5;3,6;...
%     4,1;4,2;4,3;4,4;4,5;4,6;5,1;5,2;5,3;5,4;5,5;5,6;6,1;6,2;6,3;6,4;6,5;6,6];

% dof = [1,1;1,3;1,5;2,2;2,4;2,6;3,1;3,3;3,5;4,2;4,4;4,6;5,1;5,3;5,5;6,2;6,4;6,6;];

dof = [1,1;2,2;3,3;4,4;5,5;6,6];

count = 0;
for kdof = 1:size(dof,1)
    
    count = count + 1;
    
    dof1 = dof(kdof,1);
    dof2 = dof(kdof,2);
    
    strdof(1).str = 'surge';
    strdof(2).str = 'sway';
    strdof(3).str = 'heave';
    strdof(4).str = 'roll';
    strdof(5).str = 'pitch';
    strdof(6).str = 'yaw';
    
    for k1 = 1:lw
        A11(k1) = data.hydro.A11(dof1,dof2,k1);
        A12(k1) = data.hydro.A12(dof1,dof2,k1);
        A21(k1) = data.hydro.A21(dof1,dof2,k1);
        A22(k1) = data.hydro.A22(dof1,dof2,k1);
        B11(k1) = data.hydro.B11(dof1,dof2,k1);
        B12(k1) = data.hydro.B12(dof1,dof2,k1);
        B21(k1) = data.hydro.B21(dof1,dof2,k1);
        B22(k1) = data.hydro.B22(dof1,dof2,k1);
    end
    
    figure(count)
    subplot(2,1,1)
    plot(w,A11,'-b','LineWidth',2)
    hold on
    plot(w,A12,'--k','LineWidth',2)
    plot(w,A21,'-.g','LineWidth',2)
    plot(w,A22,'-r','LineWidth',2)
    ttl = ['Coupled data - added mass (' strdof(dof1).str '-' strdof(dof2).str ')'];
    title(ttl)
    xlabel('\omega (rad/s)')
    ylabel('Added mass (kg)')
    legend('A11','A12','A21','A22')
    grid minor
    subplot(2,1,2)
    plot(w,B11,'-b','LineWidth',2)
    hold on
    plot(w,B12,'--k','LineWidth',2)
    plot(w,B21,'-.g','LineWidth',2)
    plot(w,B22,'-r','LineWidth',2)
    ttl = ['Coupled data - potential damping (' strdof(dof1).str '-' strdof(dof2).str ')'];
    title(ttl)
    xlabel('\omega (rad/s)')
    ylabel('Potential damping (kg/s)')
    legend('B11','B12','B21','B22')
    grid minor    
end

load vlcc7.txt -ascii
load suezmax17_50.txt -ascii

for k1 = 1:18
    A11_iso(19-k1) = vlcc7(18*k1,4)*1025;
    B11_iso(19-k1) = vlcc7(18*k1,5)*1025;
    A22_iso(19-k1) = suezmax17_50(18*k1,4)*1025;
    B22_iso(19-k1) = suezmax17_50(18*k1,5)*1025;
end



