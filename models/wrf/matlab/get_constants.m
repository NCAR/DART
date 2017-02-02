function [ Cp, Rd, gamma, Rv, L_c, g, T0, p0] = get_constants()
%% get_constants -  Ideally, this would take netcdf filename as input, and
% read required constants from file.  At present, just a
% repository for hardwired constants.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%% Useful constants

Rd    = 287.0;
Cp    = 7.0*Rd/2.0;
gamma = Cp / (Cp - Rd) ;
Rv    = 461; 
g     = 9.81; 
L_c   = 2.25e6; 
T0    = 300; 
p0    = 1000.e2;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
