function rh = calculate_rh(th, p, exbm, qv)
%% rh = calculate_rh(th, p, exbm, qv)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

  R  = 287;        % Gas constant for dry air
  Rv = 462;        % Gas constant for water vapor
  Cp = 1004;       % Specific heat of air (const. P)
  ep = R/Rv;       % "epsilon" - ratio of Rs
  
  tens = repmat(10,size(th));  % to easily do 10^(matrix)
  P00  = 1000;                 % Reference pressure
  
  % Calculate the pressure from the exner function
  totalp = p + exbm;
  P = P00 * (totalp) .^ (Cp/R);

  % Calculate air temperature from potential temperature - don't
  % forget to convert to Celsius
  T = th .* (P / P00) .^ (R/Cp);
  T = T - 273;

  % Calculate saturation vapor pressure - this formula comes from
  % Kerry Emanuel's "Atmospheric Convection" textbook.  T is in
  % celsius and es is in mb.
  es = 6.112 * tens .^ ((7.5 * T) ./ (T + 237.7));

  % Calculate saturation mixing ratio.  This again comes from
  % Kerry's book
  qvsat = ep * es ./ (P - es);
  
  % Relative humidity (by definition)
  rh = qv ./ qvsat;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
