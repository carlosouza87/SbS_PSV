clear variables global;close all;clc

% gen_input
% This routine should be executed in the same folder as the WAMIT data. 
% The output is a "caseid.mat" file, where "caseid" is the case identifier
% specified by the user.
% The necessary input data are:
% WAMIT output files:
% - caseid.1 files - added mass and radiation damping
% - caseid.3 files - diffraction wave loads transfer functions
% - caseid.9 files - wave drift coefficients (pressure integration method)
% - vlccT.out and suexmaxT.out files, where T is the correspondent draft 
%   for each ship.
% - MRB_vlccT.txt and MRB_suezmaxT.txt files, where T is the correspondent
%   draft for each ship. The content of the file is a 6X6 rigid body
%   matrix, e.g.:
%   7.90624E+07  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
%   0.00000E+00  7.90624E+07  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
%   0.00000E+00  0.00000E+00  7.90624E+07  0.00000E+00  0.00000E+00  0.00000E+00
%   0.00000E+00  0.00000E+00  0.00000E+00  2.59715E+10  0.00000E+00  0.00000E+00
%   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  5.06176E+11  0.00000E+00
%   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  5.06176E+11
 

%% Draft combinations for each loading case
caso(1).calado_vlcc = 7;
caso(1).calado_shuttle = 17.5;
caso(2).calado_vlcc = 14;
caso(2).calado_shuttle = 8;
caso(3).calado_vlcc = 14;
caso(3).calado_shuttle = 17.5;
caso(4).calado_vlcc = 21;
caso(4).calado_shuttle = 8;

%% File names (isolated ships)
caso(1).arq_isolado_vlcc = 'vlcc7';
caso(1).arq_isolado_shuttle = 'suezmax17_50';
caso(2).arq_isolado_vlcc = 'vlcc14';
caso(2).arq_isolado_shuttle = 'suezmax8_00';
caso(3).arq_isolado_vlcc = 'vlcc14';
caso(3).arq_isolado_shuttle = 'suezmax17_50';
caso(4).arq_isolado_vlcc = 'vlcc21';
caso(4).arq_isolado_shuttle = 'suezmax8_00';


%% Read data calculated by NYXknowledge
for k1_out = 1:4
    for k2_out = 1:length(caso(k1_out).amort)
        display(['Caso ' num2str(k1_out) ', D = ' num2str(caso(k1_out).amort(k2_out))]) % Displays the case and damping factor on the screen
        amt = caso(k1_out).amort(k2_out);
        if log10(amt) == -Inf
            amt_str = '0';
        else
            base10 = floor(log10(amt));
            mult = amt/10^base10;
            amt_str = [num2str(mult) 'E' num2str(base10)];
        end
        file_hydro_matrices = ['cs' num2str(k1_out) 'd3m']; % String with the name of the WAMIT output 
        hydro_matrices % Routine for generating added mass and potential damping coupled matrices
        %         hydro_unc % Routine for generating added mass and potential damping uncoupled matrices
        file_waveforces_data = ['cs' num2str(k1_out) 'd3m_D_' amt_str]; % String with the name of the WAMIT output with 1st order wave loads
        waveforces_data % Routine for reading 1st and 2nd order wave loads data and generating correspondent structures
        vessel1 = caso(k1_out).arq_isolado_vlcc;
        vessel2 = caso(k1_out).arq_isolado_shuttle;
        draft1 = caso(k1_out).calado_vlcc;
        draft2 = caso(k1_out).calado_shuttle;
        %         filename = ['hydrodata_' file_waveforces_data]; % String with the name of the variable with hydrodynamic data
        filename = ['cnpq_10m']
        ship_matrices % Routine for generating the structure with hydrodynamic data to be used in the simulator
        %         ship_matrices_unc
    end
end

clear all;clc