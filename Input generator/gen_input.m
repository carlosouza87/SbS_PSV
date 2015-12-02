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
% - caseid.out file is not read by the program, but the user must input
%   some data to ship_matrices which is available there.

clear all;close all;clc

% Determine type of simulation, "isimtype": 
% isimtype = 1, for Cummins equation 
% isimtype = 2, for LF + WF superposition 
isimtype = 1;  

% Determine string with the name of the WAMIT output, without any extension
caseid = ['conjunto_10']; % Water density [kg/m^3]

% Determine time step for calculation of retardation functions, if isimtype
% == 1. In this case, this will also be the time step fpr the simulation.
dt = 0.1; % Time step [s]

% Constant values adopted in WAMIT calculations
rho = 1025; % Water density [kg/m^3]
g = 9.8; % Acceleration of gravity [kg/m^3]
ULEN = 1; % WAMIT scaling factor []

% Read added mass and radiation damping, and calculate retardation functions (if isimtype == 1).
if isimtype == 1
    dt = 0.1; % Time-step for retardation functions (should match simulation time-step) [s]
end
% % Determine couplings
% cpl_mem = [1,1;1,3;1,5;2,2;2,4;2,6;3,1;3,3;3,5;4,2;4,4;4,6;5,1;5,3;5,5;6,2;6,4;6,6];
hydro_matrices

% Read 1st and 2nd order wave loads data, for isymtype == 1, and RAOs and 
% 2nd order wave loads data, for isimtype == 2. 
waveloads

% Read ships dimensions and mechanical properties
ships_data

% Consolidate the data into structures for the simulation.
data_struc

clear all;clc