function [ Cp, Rd, gamma, Rv, L_c, g, T0, p0] = get_constants()
%
% Ideally, this would take netcdf filename as input, and
% read required constants from file.  At present, just a
% repository for hardwired constants.
%

 %--Useful constants
 Cp = 1007;
 Rd = 287;  gamma = Cp / (Cp - Rd) ;
 Rv = 461; 
 g  = 9.81; 
 L_c = 2.25e6; 
 T0 = 300; 
 p0 = 1000.e2;
