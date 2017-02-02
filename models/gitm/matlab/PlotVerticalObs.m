function obs = PlotVerticalObs
% try to plot a vertical profile of observations
% 1) run MakeManyObs.m to generate a thousand independent obs.
% 2) run obs_sequence_tool to generate a single obs_seq.in
% 3) run perfect_model_obs to apply the forward observation operator
% 4) run obs_seq_to_netcdf to create something we can query
% 5) run this script to  read and plot

%
% DART $Id$
% CREDIT: Alexey Morozov

fname         = '../work/obs_epoch_001.nc';
region        = [0 360 -90 90 -Inf Inf];
CopyString    = 'truth';
QCString      = 'Quality Control';
verbose       = 1;

% nancy ... geometric or geopotential height
ObsTypeString = 'RADIOSONDE_GEOPOTENTIAL_HGT';
ObsTypeString = 'RADIOSONDE_TEMPERATURE';
ObsTypeString = 'GND_GPS_VTEC';
ObsTypeString = 'SAT_TEMPERATURE_ION';

obs   = read_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, verbose);

verts = unique(obs.Ztyp);

for i = 1:length(verts)

   figure(i); clf; orient tall;

   % find the observations on this vertical coordinate type
   inds = find(obs.Ztyp == verts(i));

   length(inds)

   X  = obs.obs(inds);
   Y  = obs.z(inds);
   qc = obs.qc(inds);

   fprintf('lon and lat are %d   %d\n',obs.lons(inds(1)), obs.lats(inds(1)))
   fprintf('min and max QC values are %d   %d\n',min(qc(:)), max(qc(:)))

   plot(X,Y,'.')
   title(sprintf('observations on coordinate type %d',verts(i)))

end

% integer, parameter :: VERTISUNDEF       = -2
% integer, parameter :: VERTISSURFACE     = -1
% integer, parameter :: VERTISLEVEL       =  1
% integer, parameter :: VERTISPRESSURE    =  2
% integer, parameter :: VERTISHEIGHT      =  3
% integer, parameter :: VERTISSCALEHEIGHT =  4

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
