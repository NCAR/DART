%% ObsTimeCheck - This is a function to explore the spatio-temporal distribution
%                 of the observations in a netCDF file created by obs_seq_to_netcdf.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

region        = [80 175 -10 60 -Inf Inf];
CopyString    = 'NCEP BUFR observation';
QCString      = 'NCEP QC index';
maxgoodQC     = 99;
verbose       = 1;
twoup         = 0;

clf;
orient landscape

for i = 1:168

   fname = sprintf('obs_epoch_%03d.nc',i);

   fignum = mod(i,8);
   if fignum == 0, fignum = 8; end

   figure(fignum)

   bob = plot_obs_netcdf(fname, 'MARINE_SFC_TEMPERATURE', region, ...
         CopyString, QCString, maxgoodQC, verbose, twoup);

   view(0,90)
   disp('Pausing, hit any key to continue ...')
   pause

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
