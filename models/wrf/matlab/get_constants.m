function [ Cp, Rd, gamma, Rv, L_c, g, T0, p0] = get_constants()
%
% Ideally, this would take netcdf filename as input, and
% read required constants from file.  At present, just a
% repository for hardwired constants.

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

 %--Useful constants
 Rd = 287.0;
 Cp = 7.0*Rd/2.0;
 gamma = Cp / (Cp - Rd) ;
 Rv = 461; 
 g  = 9.81; 
 L_c = 2.25e6; 
 T0 = 300; 
 p0 = 1000.e2;
