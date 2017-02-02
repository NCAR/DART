%% ObsPlot - plots a series of variables - over a series of times.
%            The input files are netcdf files created by obs_seq_to_netcdf.
%            The list of observations types must be hand-edited to reflect
%            the observation types of interest.
%
%            The intent is that you have a gigantic screen and can tile
%            a bunch of figure windows - one pertaining to each netcdf file.
%            If you do it right - you get a series of frames and can update
%            each frame one-at-a-time.
%
%            As it is - this is set up to create 8 figures. If each figure
%            is 3 hours ... you get one day per view.
%
%   This is not a function - the contents must be edited. The original is
%   under svn control, so you can always 'svn revert ObsPlot.m' and be good.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$


ObsTypes = {'RADIOSONDE_U_WIND_COMPONENT', ...
            'RADIOSONDE_V_WIND_COMPONENT', ...
            'RADIOSONDE_TEMPERATURE', ...
            'RADIOSONDE_SPECIFIC_HUMIDITY', ...
            'AIRCRAFT_U_WIND_COMPONENT', ...
            'AIRCRAFT_V_WIND_COMPONENT', ...
            'AIRCRAFT_TEMPERATURE', ...
            'ACARS_U_WIND_COMPONENT', ...
            'ACARS_V_WIND_COMPONENT', ...
            'ACARS_TEMPERATURE', ...
            'MARINE_SFC_U_WIND_COMPONENT', ...
            'MARINE_SFC_V_WIND_COMPONENT', ...
            'MARINE_SFC_TEMPERATURE', ...
            'MARINE_SFC_SPECIFIC_HUMIDITY', ...
            'SAT_U_WIND_COMPONENT', ...
            'SAT_V_WIND_COMPONENT', ...
            'RADIOSONDE_SURFACE_ALTIMETER', ...
            'MARINE_SFC_ALTIMETER', ...
            'LAND_SFC_ALTIMETER', ...
            'VORTEX_LAT', ...
            'VORTEX_LON', ...
            'VORTEX_PMIN', ...
            'VORTEX_WMAX', ...
            'SAT_U_WIND_COMPONENT', ...
            'SAT_V_WIND_COMPONENT'};

fname         = 'obs_epoch_003.nc';
region        = [0 360 -90 90 -Inf Inf];
region        = [80 175 -10 60 -Inf Inf];
CopyString    = 'NCEP BUFR observation';
QCString      = 'NCEP QC index';
maxgoodQC     = 99;
verbose       = 0;
twoup         = 0;

for j = 1:length(ObsTypes)
for i = 1:16     % or 168 ... or ...

   fname = sprintf('obs_epoch_%03d.nc',i);

   fignum = mod(i,8);
   if fignum == 0, fignum = 8; end

   figure(fignum)

%  bob = plot_obs_netcdf(fname, 'VORTEX_PMIN', region, ...
   bob = plot_obs_netcdf(fname, ObsTypes{j}, region, ...
         CopyString, QCString, maxgoodQC, verbose, twoup);

   view(0,90)
   fprintf('Pausing on %s, hit any key to continue ...',fname)
   pause

end
end


% ObsTypeString = 'ALL';
% for i = 1:length(ObsTypes)
%    fprintf('Observation Type %s\n', ObsTypes{i})
%    bob = plot_obs_netcdf(fname, ObsTypes{i}, region, ...
%          CopyString, QCString, maxgoodQC, verbose, twoup);
%    view(0,90)
%    disp('Pausing, hit any key to continue ...')
%    pause
% end


% 1  'RADIOSONDE_U_WIND_COMPONENT'
% 2  'RADIOSONDE_V_WIND_COMPONENT'
% 3  'RADIOSONDE_TEMPERATURE'
% 4  'RADIOSONDE_SPECIFIC_HUMIDITY'
% 5  'AIRCRAFT_U_WIND_COMPONENT'
% 6  'AIRCRAFT_V_WIND_COMPONENT'
% 7  'AIRCRAFT_TEMPERATURE'
% 8  'ACARS_U_WIND_COMPONENT'
% 9  'ACARS_V_WIND_COMPONENT'
% 10 'ACARS_TEMPERATURE'
% 11 'MARINE_SFC_U_WIND_COMPONENT'
% 12 'MARINE_SFC_V_WIND_COMPONENT'
% 13 'MARINE_SFC_TEMPERATURE'
% 14 'MARINE_SFC_SPECIFIC_HUMIDITY'
% 15 'SAT_U_WIND_COMPONENT'
% 16 'SAT_V_WIND_COMPONENT'
% 17 'RADIOSONDE_SURFACE_ALTIMETER'
% 18 'MARINE_SFC_ALTIMETER'
% 19 'LAND_SFC_ALTIMETER'
% 20 'VORTEX_LAT'
% 21 'VORTEX_LON'
% 22 'VORTEX_PMIN'
% 23 'VORTEX_WMAX'
% 24 'SAT_U_WIND_COMPONENT'
% 25 'SAT_V_WIND_COMPONENT'

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
