% clear all;close all;clc
clear variables global;close all;clc

betaw_list(1,1) = 195;
betaw_list(2,1) = 315;
Hs_list(1,1) = 1.5;
Hs_list(2,1) = 2.5;
Hs_list(3,1) = 3.5;
Tp_list(1,1) = 8.0;
Tp_list(2,1) = 11.5;
Tp_list(3,1) = 15.0;
strmoor = ['moor846';'moor844';'moor644'];

count = 0;
for moor_arr = 1:3
    for caso = 1:3
        caso
        for k_incid = 1:2
            betaw = betaw_list(k_incid);
            for k_Hs = 1:3
                Hs = Hs_list(k_Hs);
                for k_Tp = 1:3
                    Tp = Tp_list(k_Tp);
                    % wind
                    gammaw = betaw;
                    Uw = 7.7; % wind velocity [m/s]
                    % current
                    alphac = pi;    % current incidence direction [deg]
                    Uc= 0;    % current velocity [m/s]
                    main
                    simname = ['sim_caso' num2str(caso) '_Hs' num2str(Hs*10) 'Tp' num2str(Tp*10) 'incid' num2str(betaw) strmoor(moor_arr,:)];
                    save(simname,'tsim','data','variable')
                    clear variable data y
                    count = count + 1
                end
            end
        end
    end
end


