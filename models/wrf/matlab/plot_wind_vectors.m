function obs = plot_wind_vectors(fname,ncname,platform)
% plot_wind_vectors(fname,ncname,platform)
%
%  102 RADIOSONDE
% 2223 SAT
% 4243 BUOY
% 4748 SHIP
% 5758 AIREP
% 6768 PILOT
% 8182 QKSWND
%
% fname = 'wind_vectors.006.dat';
% ncname = 'obs_diag_output.nc';
% platform = 8182;    % QKSWND
% ilevel = 1;         % surface == 1?
%
% plot_wind_vectors(fname,ncname,platform)

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

f = netcdf(ncname,'nowrite');
platforms = f{'ObservationTypes'}(:);
timeunits = f{'time'}.units(:);
close(f);

timebase     = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin   = datenum(timebase(1),timebase(2),timebase(3));

obs    = Load_Parse(fname, platforms, platform );

obs.times = obs.times + timeorigin;
t1 = datestr(min(obs.times));
t2 = datestr(max(obs.times));

clf;
axlims = DrawBackground(obs);
h1     = quiver(obs.lon, obs.lat, obs.U, obs.V);
h2     = plot(obs.lon, obs.lat, '.');
set(h2,'MarkerSize',8)
hold off;

h = title({sprintf('%s %s %s',t1,platform,t2),'all levels - including 99s bad UVs'});
set(h,'FontSize',18)

xlabel(fname)


function axlims = DrawBackground(obs, platform)
%----------------------------------------------------------------------

% Figure out bounds of the data
axlims = [ min(obs.lon) max(obs.lon) min(obs.lat) max(obs.lat) ];

% It is nice to have a little padding around the perimeter
dx = 0.025 * (axlims(2) - axlims(1));
dy = 0.025 * (axlims(4) - axlims(3));

axis(axlims + [-dx dx -dy dy])

% It is nice to know where the land is
worldmap('light');
hold on;




function obs = Load_Parse(fname, platforms, platform)
%----------------------------------------------------------------------
% Makes no attempt to find/replace/identify MISSING values

obsmat = load(fname);

platformIDs = unique(obsmat(:,1));
uid         = floor(platformIDs/100);
vid         = platformIDs - uid*100;
Ustrings    = platforms(uid,:);
Vstrings    = platforms(vid,:);

nplatforms = length(uid);
pid        = [];

for i = 1:nplatforms
   uindex = findstr(Ustrings(i,:),'_U_WIND_COMPONENT') - 1;
   vindex = findstr(Vstrings(i,:),'_V_WIND_COMPONENT') - 1;

   Ubase  = Ustrings(i,1:uindex);
   Vbase  = Vstrings(i,1:vindex);

   if ~strcmp(Ubase,Vbase)
      error('U and V wind component platforms do not match')
   end

   disp(sprintf('Available platforms in %s are %s',fname,Ubase))

   if strcmp(lower(Ubase),lower(platform))
      pid = platformIDs(i);
   end

end

if isempty(pid)
   error('no %s observations in %s',platform,fname)
end

inds   = find(obsmat(:,1) == pid);

if isempty(inds)
   error('no observations of type %d in %s',pid,fname)
end

obs.platform = obsmat(inds, 1);
obs.day      = obsmat(inds, 2);
obs.seconds  = obsmat(inds, 3);
obs.lon      = obsmat(inds, 4);
obs.lat      = obsmat(inds, 5);
obs.level    = obsmat(inds, 6);
obs.levind   = obsmat(inds, 7);
obs.Uqc      = obsmat(inds, 8);
obs.Vqc      = obsmat(inds, 9);
obs.U        = obsmat(inds,10);
obs.V        = obsmat(inds,11);
obs.Upr      = obsmat(inds,12);
obs.Vpr      = obsmat(inds,13);
obs.Upo      = obsmat(inds,14);
obs.Vpo      = obsmat(inds,15);
obs.times = obs.day + obs.seconds/86400;

clear obsmat
