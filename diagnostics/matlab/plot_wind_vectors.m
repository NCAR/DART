function data = plot_wind_vectors( fname, platform, CopyString, QCString, varargin )
%% plot_wind_vectors creates maps with wind vector images overlain.
%
% The required arguments are:
% fname    ... the output netCDF file from obs_seq_to_netcdf
% platform ... a string to represent the observation platform.
%              usually 'RADIOSONDE', 'SAT', 'METAR', ... get hints from:
%              ncdump -v ObsTypesMetaData obs_epoch_xxx.nc | grep _U_  
% CopyString ... which observation copy is of interest? 
%              ncdump -v CopyMetaData obs_epoch_xxx.nc
% QCString   ... which QC copy is of interest?
%              ncdump -v QCMetaData obs_epoch_xxx.nc
% 
% The optional arguments are:
% region      ... specifies that horizontal & vertical area of interest. If not present,
%                 all available observations are used. 
% scalefactor ... provides control over the plotted size of the 
%                 wind vectors. A smaller number results in a 
%                 bigger wind vector. If not present, a value of 10.0 is
%                 used.
%
% EXAMPLE 1:
% fname        = 'obs_epoch_001.nc';
% platform     = 'SAT';    % usually 'RADIOSONDE', 'SAT', 'METAR', ...
% CopyString   = 'NCEP BUFR observation';
% QCString     = 'DART quality control';
% region       = [0 360 0 90 1020 500];    % 
% scalefactor = 5;     % reference arrow magnitude
%
% obs = plot_wind_vectors(fname, platform, CopyString, QCString, ...
%         'region', region, 'scalefactor', scalefactor);
%
%
% EXAMPLE 2 (CONUS domain):
%
% region = [210 310 12 65 -Inf Inf];
% obs = plot_wind_vectors('obs_epoch_001.nc', 'SAT', ...
%       'NCEP BUFR observation', 'DART quality control','region',region);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Set sensible defaults

region      = [0 360 -90 90 -Inf Inf];
scalefactor = 10;

% Harvest input

if nargin ~= 4 
   [region, scalefactor] = parseinput(varargin{:});
end

%% Start the ball rolling

if (exist(fname,'file') ~= 2) 
   error('%s does not exist',fname)
end

data.filename    = fname;
data.platform    = platform;
data.copystring  = CopyString;
data.qcstring    = QCString;
data.region      = region;
data.scalefactor = scalefactor;

%% Read the observation sequence

[UtypeString, VtypeString] = FindObsType( fname, platform );
verbose = 0;
Uobs = read_obs_netcdf(fname, UtypeString, region, CopyString, QCString, verbose);
Vobs = read_obs_netcdf(fname, VtypeString, region, CopyString, QCString, verbose);

if (length(Uobs.obs) ~= length(Vobs.obs))
   error('Houston, we have a problem.')
end

lonmismatch = Uobs.lons ~= Vobs.lons;
latmismatch = Uobs.lats ~= Vobs.lats;
  zmismatch = Uobs.z    ~= Vobs.z;

if ( sum(lonmismatch) ~= 0),
   warning('DART:UVcollocaton','There are %d mismatched (in longitude) observations',sum(lonmismatch))
end 
if ( sum(latmismatch) ~= 0),
   warning('DART:UVcollocaton','There are %d mismatched (in latitude) observations',sum(latmismatch))
end 
if ( sum(zmismatch) ~= 0),
   warning('DART:UVcollocaton','There are %d mismatched (in vertical) observations',sum(zmismatch))
end 

if (sum(lonmismatch) == length(lonmismatch))
   clf;
   axis([0 1 0 1]); axis image
   h = text(0.5,0.5,sprintf('%s has no %s data',data.filename,data.platform));
   set(h,'Interpreter','none')
   return
end

%% must only use the observations that are co-located to 
%  set up the plotting structure.

inds = ((Uobs.lons == Vobs.lons) & ...
        (Uobs.lats == Vobs.lats) & ...
        (Uobs.z    == Vobs.z));

data.time = Uobs.time(inds);
data.lon  = Uobs.lons(inds);
data.lat  = Uobs.lats(inds);
data.z    = Uobs.z(inds);
data.Uqc  = Uobs.qc(inds);
data.Vqc  = Vobs.qc(inds);
data.U    = Uobs.obs(inds);
data.V    = Vobs.obs(inds);
data.level1    = min(Uobs.z);
data.levelN    = max(Uobs.z);
data.levstring = sprintf('%.2f to %.2f',data.level1,data.levelN);

%% Start the Plotting

clf;
axish = gca;
axlims = DrawBackground( data );
set(axish,'Layer','top')

if ( isfinite(strfind(lower(QCString),'dart')) )
   % We know how to interpret QC codes
   goodUV = find( (data.Uqc < 2) & (data.Vqc < 2));
   baadUV = find( (data.Uqc > 1) & (data.Vqc > 1));
   goodU  = find( (data.Uqc < 2) & (data.Vqc > 1));
   goodV  = find( (data.Uqc > 1) & (data.Vqc < 2));
else
   % We do not know how to interpret QC codes, so they
   % are all 'good'
   baadUV = [];
   goodU  = [];
   goodV  = [];
end

legh   = [];
legstr = {};

if ~ isempty(goodUV)
   hgood  = obs2plot(data, goodUV, [0 0 0] );
   legh   = [legh; hgood];
   legstr{length(legstr)+1} = sprintf('%d ''good''',length(goodUV));
end

if ~ isempty(baadUV)
   hbaadUV = obs2plot(data, baadUV, [1 0 0] );
   legh    = [legh; hbaadUV];
   legstr{length(legstr)+1} = sprintf('%d ''badU  badV''',length(baadUV));
end

if ~ isempty(goodU)
   hgoodU = obs2plot(data, goodU, [0 1 0] );
   legh   = [legh; hgoodU];
   legstr{length(legstr)+1} = sprintf('%d ''goodU badV''',length(goodU));
end

if ~ isempty(goodV)
   hgoodV = obs2plot(data, goodV, [0 0 1] );
   legh   = [legh; hgoodV];
   legstr{length(legstr)+1} = sprintf('%d ''badU  goodV''',length(goodV));
end

t1 = datestr(min(data.time),'yyyy-mm-dd HH:MM:SS');
t2 = datestr(max(data.time),'yyyy-mm-dd HH:MM:SS');
h = title({sprintf('%s %s %s',t1,platform,t2), ...
           sprintf('levels %s ',data.levstring)});
set(h,'FontSize',18)

h = xlabel(data.filename); set(h,'Interpreter','none');

legend(legh,legstr,'Location','Best','FontSize',18)

hold off;


function axlims = DrawBackground( obs )
%======================================================================

if  isempty(obs.region)
   % Figure out bounds of the data
   axlims = [ min(obs.lon) max(obs.lon) min(obs.lat) max(obs.lat) ] ;
else
   axlims = obs.region(1:4);
end

% It is nice to have a little padding around the perimeter
dx = 0.05 * (axlims(2) - axlims(1));
dy = 0.05 * (axlims(4) - axlims(3));
axlims(1:4) = axlims(1:4) + [-dx dx -dy dy];
axis(axlims)

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



function [ustring, vstring] = FindObsType( ncname, platform )
%======================================================================
% Makes no attempt to find/replace/identify MISSING values
%
% data.filename     incoming filename with the data
% data.ncname       incoming filename with the metadata
% data.platform     the observation platform of interest
% data.region       the region of interest [xmin xmax ymin ymax]
% data.scalefactor  the reference wind vector magnitude
% data.platforms    the observation platforms in the incoming file
% data.timeunits    the units string e.g. 'days since yyyy-mm-dd'
% data.timeorigin

ObsTypeStrings = ncread(ncname,'ObsTypesMetaData')';

% Find the types of data in the obs_sequence file for this epoch.
% Turns out all the winds are either xxxx_U_WIND_COMPONENT or xxxx_U_10_METER_WIND
% so either way, they start out xxxx_U_

utarget = sprintf('%s_U_',strtrim(platform));
vtarget = sprintf('%s_V_',strtrim(platform));
n       = length(utarget);
ustring = [];
vstring = [];

for i = 1:size(ObsTypeStrings,1)

   utf = strncmpi(ObsTypeStrings(i,:),utarget,n);
   vtf = strncmpi(ObsTypeStrings(i,:),vtarget,n);

   if ( utf )
      ustring = deblank(ObsTypeStrings(i,:));
      uindex  = i;
   end 
   if ( vtf )
      vstring = deblank(ObsTypeStrings(i,:));
      vindex  = i;
   end 

end

if ( isempty(ustring) || isempty(vstring) )
   error('no %s winds in %s',platform,ncname)
end

% echo a little informational statement about the number of obs

obs_type = ncread(ncname,'obs_type');
numU     = sum(obs_type == uindex);
numV     = sum(obs_type == vindex);
fprintf('%8d %s observations in %s\n', numU, ustring, ncname)
fprintf('%8d %s observations in %s\n', numV, vstring, ncname)

if (numU ~= numV)
   error('Different number of U,V observations. Dying ...')
end



function [region, scalefactor] = parseinput(varargin)
%======================================================================
% try to parse the input pairs ... which must be pairs

if (mod(length(varargin),2) ~= 0)
   error('Wrong number (%d) of optional arguments. Must be parameter/value pairs: ''region'',[0 360 -90 90 1020 500]',length(varargin)) 
end

npairs = length(varargin)/2;

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

if (length(region) < 4)
   region = [0 360 -90 90 -Inf Inf];
elseif (length(region) == 4)
   region = [region -Inf Inf];
elseif (length(region) ~= 6)
   warning('DART:region input','region must be length 4 or 6') 
   error('Unable to interpret region - %s',num2str(region))
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
