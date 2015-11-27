clear variables global;close all;clc

% gen_input
% This routine should be executed in the same folder as the WAMIT data.
% The output is a "caseid.mat" file, where "caseid" is the case identifier
% specified by the user.
% The necessary input data are:
% WAMIT output files:
% - caseid.1 files - added mass and radiation damping;
% - caseid.3 files - diffraction wave loads transfer functions (isimtype=1), OR
% - caseid.4 files - response amplitude operators (isimtype = 2);
% - caseid.9 files - wave drift coefficients (pressure integration method);
% - vlccT.out and suexmaxT.out files, where T is the correspondent draft for each ship;
% - MRB_vlccT.txt and MRB_suezmaxT.txt files, where T is the correspondent
%   draft for each ship. The content of the file is a 6X6 rigid body
%   matrix, e.g.:
%   7.90624E+07  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
%   0.00000E+00  7.90624E+07  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
%   0.00000E+00  0.00000E+00  7.90624E+07  0.00000E+00  0.00000E+00  0.00000E+00
%   0.00000E+00  0.00000E+00  0.00000E+00  2.59715E+10  0.00000E+00  0.00000E+00
%   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  5.06176E+11  0.00000E+00
%   0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  5.06176E+11


% Determine type of simulation, "isimtype": 
% isimtype = 1, for Cummins equation 
% isimtype = 2, for LF + WF superposition 
isimtype = 1;  

% Determine string with the name of the WAMIT output, without any extension
caseid = ['conjunto']; % Water density [kg/m^3]

% Constant values adopted in WAMIT calculations
rho = 1025; % Water density [kg/m^3]
g = 9.8; % Acceleration of gravity [kg/m^3]
ULEN = 1; % WAMIT scaling factor []

% Read added mass and radiation damping, and calculate retardation functions (if isimtype == 1).
if isimtype == 1
    dt = 0.1; % Time-step for retardation functions (should match simulation time-step) [s]
end
hydro_matrices

% Read 1st and 2nd order wave loads data, for isymtype == 1, and RAOs and 
% 2nd order wave loads data, for isimtype == 2. 
waveforces_data

% Define rigid-body matrices, read further information from WAMIT .out files and
% consolidate the data into structures for the simulation.
ship_matrices

clear all;clc