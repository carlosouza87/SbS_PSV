function [mu,chip] = ssmemory(Ar,Br,Cr,chi,u)
% ssmemory(Ar,Br,Cr,Dr,chi,dnu)
% Performs state-space approximation of the retardation function
% 
% Ar, Br and Cr are the state-space matrixes
% chi is the state vector
% u is the input
% 
% Outputs are:
% mu is the retardation function approximation
% chip is the time derivative of the state vector

chip = Ar*chi + Br*u;
mu = Cr*chi;