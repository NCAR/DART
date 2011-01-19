function obs = plot_wind_vectors( fname, ncname, platform, varargin )
%% plot_wind_vectors creates maps with wind vector images overlain.
%
% The required arguments are:
% fname    ... the file name containing the wind observation pairs
% ncname   ... the output netCDF file from obs_diag.f90 (for some metadata)
% platform ... a string to represent the observation platform.
% 
% The optional arguments are:
% levels      ... specifies the vertical area of interest. If not present,
%                 all vertical levels are used.
% region      ... specifies that horizontal area of interest. If not present,
%                 all available observations are used. 
% scalefactor ... provides control over the plotted size of the 
%                 wind vectors. A smaller number results in a 
%                 bigger wind vector. If not present, a value of 10.0 is
%                 used.
%
% fname    = 'wind_vectors.006.dat';
% ncname   = 'obs_to_table_output.nc';
% platform = 'RADIOSONDE';
% levels   = [1020 500];
% region   = [0 360 0 90];    % 
% scalefactor = 5;     % reference arrow magnitude
%
% obs = plot_wind_vectors(fname, ncname, platform, ...
%         'levels', levels, 'region', region, 'scalefactor', scalefactor);
%

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

narg = nargin;

if narg == 3 

   levels   = [];
   region   = [];
   scalefactor = 10;

else

   [levels, region, scalefactor] = parseinput(varargin{:});

end

data.filename    = fname;
data.ncname      = ncname;
data.platform    = platform;
data.levels      = levels;
data.region      = region;
data.scalefactor = scalefactor;

f = netcdf(ncname,'nowrite');
data.platforms = f{'ObservationTypes'}(:);
data.timeunits = f{'time_bounds'}.units(:);
close(f);

timebase   = sscanf(data.timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
data.timeorigin = datenum(timebase(1),timebase(2),timebase(3));
obs        = Load_Parse(data);

if  isempty(obs)
   clf;
   axis([0 1 0 1]); axis image
   h = text(0.5,0.5,sprintf('%s has no %s data',data.filename,data.platform));
   set(h,'Interpreter','none')
   return
end

obs.times  = obs.times + data.timeorigin;
t1         = datestr(min(obs.times));
t2         = datestr(max(obs.times));

clf;
axlims     = DrawBackground( obs );


goodUV = find( (obs.Uqc < 1) & (obs.Vqc < 1));
baadUV = find( (obs.Uqc > 1) & (obs.Vqc > 1));
goodU  = find( (obs.Uqc < 1) & (obs.Vqc > 0));
goodV  = find( (obs.Vqc < 1) & (obs.Uqc > 0));

legh   = [];
legstr = {};

if ~ isempty(goodUV)
   hgood  = obs2plot(obs, goodUV, [0 0 0] );
   legh   = [legh; hgood];
   legstr{length(legstr)+1} = sprintf('%d good',length(goodUV));
end

if ~ isempty(baadUV)
   hbaadUV = obs2plot(obs, baadUV, [1 0 0] );
   legh    = [legh; hbaadUV];
   legstr{length(legstr)+1} = sprintf('%d badUbadV',length(baadUV));
end

if ~ isempty(goodU)
   hgoodU = obs2plot(obs, goodU, [0 1 0] );
   legh   = [legh; hgoodU];
   legstr{length(legstr)+1} = sprintf('%d goodUbadV',length(goodU));
end

if ~ isempty(goodV)
   hgoodV = obs2plot(obs, goodV, [0 0 1] );
   legh   = [legh; hgoodV];
   legstr{length(legstr)+1} = sprintf('%d badUgoodV',length(goodV));
end

h = title({sprintf('%s %s %s',t1,platform,t2), ...
           sprintf('levels %s ',obs.levelstring)});
set(h,'FontSize',18)

h = xlabel(data.filename); set(h,'Interpreter','none');

legend(legh,legstr,'Location','NorthWestOutside')

hold off;



function axlims = DrawBackground( obs )
%======================================================================

if  isempty(obs.region)
   % Figure out bounds of the data
   axlims = [ min(obs.lon) max(obs.lon) min(obs.lat) max(obs.lat) ] ;
else
   axlims = obs.region;
end

% It is nice to have a little padding around the perimeter
dx = 0.05 * (axlims(2) - axlims(1));
dy = 0.05 * (axlims(4) - axlims(3));

axis(axlims + [-dx dx -dy dy])
axis image

% It is nice to know where the land is
worldmap('light');
hold on;

% Plot a wind vector pair of known length for reference
tx = axlims(1)+dx;
ty = axlims(4)-dy;

U = 10/obs.scalefactor; 
V =  0/obs.scalefactor;
h = quiver(tx, ty, U, V, 0.0 ,'LineWidth',4.0,'Color','k');
set(h,'LineWidth',3.0)

U =  0/obs.scalefactor; V = 10/obs.scalefactor;
h = quiver(tx, ty, U, V, 0.0 ,'LineWidth',4.0,'Color','k');
set(h,'LineWidth',3.0)
h = text(tx, ty-0.1*dy,'10 m/s','VerticalAlignment','top');



function h1 = obs2plot(obs, mask, colspec )
%======================================================================
% Actually plots the wind vectors. to defeat Matlab's automatic 
% scaling, we use a scalefactor of 0.0 and manually scale the
% data from the user input (or default).

lon    = obs.lon(mask);
lat    = obs.lat(mask);
U      = obs.U(mask)/obs.scalefactor;
V      = obs.V(mask)/obs.scalefactor;
h1     = quiver(lon, lat, U, V, 0.0);
h2     = plot(lon, lat, '.','MarkerSize',4);
set(h1,'Color',colspec)
set(h2,'Color',colspec)



function obs = Load_Parse( data )
%======================================================================
% Makes no attempt to find/replace/identify MISSING values
%
% data.filename     incoming filename with the data
% data.ncname       incoming filename with the metadata
% data.platform     the observation platform of interest
% data.levels       the top/bottom levels of interest
% data.region       the region of interest [xmin xmax ymin ymax]
% data.scalefactor  the reference wind vector magnitude
% data.platforms    the observation platforms in the incoming file
% data.timeunits    the units string e.g. 'days since yyyy-mm-dd'
% data.timeorigin

obsmat = load(data.filename);

% Find the types of data in the obs_sequence file for this epoch.

platformIDs = unique(obsmat(:,1));
uid         = floor(platformIDs/100);
vid         = platformIDs - uid*100;
Ustrings    = data.platforms(uid,:);
Vstrings    = data.platforms(vid,:);

nplatforms = length(uid);
pid        = [];
obs        = [];
levelstring = [];
regionstring = [];

% This block divines the platform string and companion obs_kinds 
% from the netCDF file metadata - the only reason we need the netCDF file.

for i = 1:nplatforms
   uindex = findstr(Ustrings(i,:),'_U_WIND_COMPONENT') - 1;
   vindex = findstr(Vstrings(i,:),'_V_WIND_COMPONENT') - 1;

   if (isempty(uindex) | isempty(vindex)) 
      uindex = findstr(Ustrings(i,:),'_U_10_METER_WIND') - 1;
      vindex = findstr(Vstrings(i,:),'_V_10_METER_WIND') - 1;
   end 

   Ubase  = Ustrings(i,1:uindex);
   Vbase  = Vstrings(i,1:vindex);

   if ~strcmp(Ubase,Vbase)
      error('U and V wind component platforms do not match')
   end

   % Determine what numeric pid corresponds to the desired platform string
   if strcmp(lower(Ubase),lower(data.platform))
      pid = platformIDs(i);
   end

   % echo a little informational statement about how many obs of
   % this type there are in this assimilation period.
   inds = find(obsmat(:,1) == platformIDs(i));
   nobs = length(inds);

   disp(sprintf('%6d %14s observations in %s (%4d)',nobs,Ubase,data.filename,platformIDs(i)))

end

% This block extracts just the desired observations based on platform.

if isempty(pid)
   disp(sprintf('no %s observations in %s', data.platform, data.filename))
   return
end

inds   = find(obsmat(:,1) == pid);

if isempty(inds)
   disp(sprintf('no %s observations (type %d) in %s',data.platform, pid, data.filename))
   return
end

platform = obsmat(inds, 1);
day      = obsmat(inds, 2);
seconds  = obsmat(inds, 3);
lon      = obsmat(inds, 4);
lat      = obsmat(inds, 5);
level    = obsmat(inds, 6);
Uqc      = obsmat(inds, 7);
Vqc      = obsmat(inds, 8);
U        = obsmat(inds, 9);
V        = obsmat(inds,10);
if ( size(obsmat,2) > 10 )
   Upr   = obsmat(inds,11);
   Vpr   = obsmat(inds,12);
   Upo   = obsmat(inds,13);
   Vpo   = obsmat(inds,14);
end
times    = day + seconds/86400;

%--------------------------------------------------
% Subset the levels of interest
%--------------------------------------------------

if ( isempty(data.levels) ) 
   % then we want all levels, do nothing ...
   level1 = min(level);
   levelN = max(level);

   levelstring = sprintf('all (%.2f to %.2f)',level1,levelN);
else
   level1 = min(data.levels);
   levelN = max(data.levels);

   levelstring = sprintf('%.2f to %.2f',level1,levelN);

   inds = find ((level >= level1) & (level <= levelN));

   if (length(inds) == 0)
      disp(sprintf('no %s observations in %s', data.platform, levelstring))
      return
   end

   platform = platform(inds);
   day      =      day(inds);
   seconds  =  seconds(inds);
   lon      =      lon(inds);
   lat      =      lat(inds);
   level    =    level(inds);
   Uqc      =      Uqc(inds);
   Vqc      =      Vqc(inds);
   U        =        U(inds);
   V        =        V(inds);
   if ( size(obsmat,2) > 10 )
      Upr   =      Upr(inds);
      Vpr   =      Vpr(inds);
      Upo   =      Upo(inds);
      Vpo   =      Vpo(inds);
   end
   times    =    times(inds);
end

%--------------------------------------------------
% Subset the region of interest
% for the moment, we are not supporting wrapping at Greenwich.
%--------------------------------------------------

if isempty(data.region)
   % then we want the entire dataset, do nothing ... 
   regionstring = 'global';
else
   lon1 = data.region(1);
   lonN = data.region(2);
   lat1 = data.region(3);
   latN = data.region(4);

   regionstring = sprintf('(%.2f -> %.2f, %.2f -> %.2f)',lon1,lonN,lat1,latN);

   inds = find ((lon >= lon1) & (lon <= lonN) & ...
                (lat >= lat1) & (lat <= latN));

   if (length(inds) == 0)
      disp(sprintf('no %s observations in %s', data.platform, regionstring))
      return
   end

   platform = platform(inds);
   day      =      day(inds);
   seconds  =  seconds(inds);
   lon      =      lon(inds);
   lat      =      lat(inds);
   level    =    level(inds);
   Uqc      =      Uqc(inds);
   Vqc      =      Vqc(inds);
   U        =        U(inds);
   V        =        V(inds);
   if ( size(obsmat,2) > 10 )
      Upr   =      Upr(inds);
      Vpr   =      Vpr(inds);
      Upo   =      Upo(inds);
      Vpo   =      Vpo(inds);
   end
   times    =    times(inds);
end

%--------------------------------------------------
% Insert NaN's for missing values
%--------------------------------------------------

U(   U   < -900) = NaN;
V(   V   < -900) = NaN;
if ( size(obsmat,2) > 10 )
   Upr( Upr < -900) = NaN;
   Vpr( Vpr < -900) = NaN;
   Upo( Upo < -900) = NaN;
   Vpo( Vpo < -900) = NaN;
end

%--------------------------------------------------
% Create the output structure.
%--------------------------------------------------

obs.platform    = platform;
obs.day         = day;
obs.seconds     = seconds;
obs.lon         = lon;
obs.lat         = lat;
obs.level       = level;
obs.Uqc         = Uqc;
obs.Vqc         = Vqc;
obs.U           = U;
obs.V           = V;
if ( size(obsmat,2) > 10 )
   obs.Upr      = Upr;
   obs.Vpr      = Vpr;
   obs.Upo      = Upo;
   obs.Vpo      = Vpo;
end
obs.times       = times;
obs.levels      = data.levels;
obs.region      = data.region;
obs.scalefactor = data.scalefactor;
obs.levelstring = levelstring;


function [levels, region, scalefactor] = parseinput(varargin)
%======================================================================
% try to parse the input pairs ... which must be pairs

if (mod(length(varargin),2) ~= 0)
   error('Wrong number (%d) of optional arguments. Must be parameter/value pairs: ''levels'',[1000 500]',length(varargin)) 
end

npairs = length(varargin)/2;

levels      = [];
region      = [];
scalefactor = 10.0;

for i = 1:2:length(varargin)
   switch lower(varargin{i})
      case 'levels'
         levels = varargin{i+1};
      case 'region'
         region = varargin{i+1};
      case 'scalefactor'
         scalefactor = varargin{i+1};
      otherwise
         error('Unknown parameter %s',varargin{i})
   end
end

% Make sure the levels array has a top/bottom
if  ~ isempty(levels)
end

% Make sure the geographic array makes sense
if  ~ isempty(region)
end

% Make sure the scalefactor makes sense
if  ~ isempty(scalefactor)
end

