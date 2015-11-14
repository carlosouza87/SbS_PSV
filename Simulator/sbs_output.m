clear global variables;close all;clc

% caso = input('Caso: ');
% moor_arr = input('Mooring arrangement: ');
betaw_list(1,1) = 195;
betaw_list(2,1) = 315;
Hs_list(1,1) = 1.5;
Hs_list(2,1) = 2.5;
Hs_list(3,1) = 3.5;
Tp_list(1,1) = 8.0;
Tp_list(2,1) = 11.5;
Tp_list(3,1) = 15.0;
strmoor = ['846';'844';'644'];

count = 0;
for caso = 1:3
    for moor_arr = 1:3
        for k_incid = 1:2
            betaw = betaw_list(k_incid);
            for k_Hs = 1:3
                Hs = Hs_list(k_Hs);
                for k_Tp = 1:3
                    Tp = Tp_list(k_Tp);
                    filename = ['Caso ' num2str(caso) '\sim_caso' num2str(caso) '_Hs' num2str(Hs*10) 'Tp' num2str(Tp*10) 'incid' num2str(betaw) 'moor' strmoor(moor_arr,:)];
                    load(filename)
                    mbl = zeros(4,1);
                    swl = zeros(4,1);
                    mean_tract = zeros(4,1);
                    std_tract = zeros(4,1);
                    if moor_arr == 1
                        kpos = [1;9;11;13];
                    elseif moor_arr == 2
                        kpos = [1;9;11;13];
                    elseif  moor_arr == 3
                        kpos = [1;7;9;11];
                    end                        
                    for k1 = 1:4
                        T = variable.sgm(kpos(k1),:)*data.mooring.D^2*pi/4;
                        mbl(k1) = sum(T>data.mooring.mbl(k1));
                        swl(k1) = sum(T>data.mooring.swl(k1));
                        mean_tract(k1) = mean(T);
                        std_tract(k1) = std(T);
                    end
                    ident = 1;
                    if sum(swl) > 0
                        ident = 2;
                    end
                    if sum(mbl) > 0
                        ident = 3;
                    end
                    if variable.break == 1
                        ident = 4;
                    end
                    mooring_results.case(caso).mooring(moor_arr).incid(k_incid).table(k_Hs,k_Tp) = ident;
                    mooring_results.case(caso).mooring(moor_arr).incid(k_incid).mean_tract{k_Hs,k_Tp} = mean_tract;
                    mooring_results.case(caso).mooring(moor_arr).incid(k_incid).std_tract{k_Hs,k_Tp} = std_tract;
                    mooring_results.case(caso).mooring(moor_arr).incid(k_incid).betaw = betaw_list(k_incid);
                    mooring_results.case(caso).mooring(moor_arr).incid(k_incid).moor_arr = strmoor(moor_arr,:);
                    mooring_results.case(caso).mooring(moor_arr).incid(k_incid).condition = caso;
                    count = count + 1
                end
            end
        end
    end
end

save analysis mooring_results