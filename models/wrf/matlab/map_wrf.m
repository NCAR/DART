function h = map_wrf(fname, varname, levelindx, timeindx, copystring )
%% map_wrf creates an image of a WRF field without using the mapping toolbox.
%
%% Example using a DART-style diagnostic file, i.e.:
%
% fname      = '/glade/proj3/image/hliu/200812new/cwb_icbc/12_01/Prior_Diag.nc';
% varname    = 'U_d01';
% copystring = 'ensemble mean';
% levelindx  = 10;
% timeindx   = 1;
% map_wrf(fname, varname, levelindx, timeindx, copystring);
% worldmap;
% axis off;
%
%
%% Example using a WRF netCDF file:
%
% fname      = '/glade/proj3/image/hliu/ICBC_from_cwb/wrfinput_d01';
% levelindx  = 10;
% timeindx   = 1;
% map_wrf(fname, 'U', levelindx, timeindx);
% worldmap;   % superimpose some low-res coastlines 
% axis off;
%
%% layer on the locations of some observations:
%
% obsfile = 'obs_epoch_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'observation';
% QCString      = 'DART quality control';
% verbose       = 1;   % anything > 0 == 'true'
% obs = read_obs_netcdf(obsfile, ObsTypeString, region, CopyString, QCString, verbose);
% hold on;
% plot(obs.lons, obs.lats, 'kd', 'MarkerFaceColor','k')

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( nargin == 3 )
   levelindx = 1;
   timeindx  = 1;
elseif (nargin == 4 ) 
   timeindx  = 1;
elseif (nargin ~= 5 ) 
   error('Wrong number of arguments ... must have 3,4, or 5')
end

% Check to ensure the file and the desired variable exist

if (exist(fname,'file') ~= 2) 
   error('%s does not exist',fname)
end

varexist(fname,{varname})

%% Get the dimension information.

xy          = GetWRFlatlon( fname, varname);
levels      = GetWRFlevels( fname, varname);
times       = GetWRFtimes(  fname, varname);
timestrings = datestr(times,31);

%% Get the variable units, description, etc.
% Determine variable to plot ...

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Attribute)
   attname = varinfo.Attribute(i).Name;
   attvalu = varinfo.Attribute(i).Value;
   switch lower(attname)
      case 'name'
         varinfo.name = attvalu;
      case 'long_name'
         varinfo.long_name = attvalu;
      case 'units'
         varinfo.units = attvalu;
   end
end

copydim = find(strncmp('copy',varinfo.Dimension,length('copy')));

if (isempty(copydim))
   copyindex = NaN;
else
   copyindex = get_copy_index(fname,copystring);
end

for itime = timeindx

   plot_title = {fname, ...
         sprintf('levelindex %d %s %s',levelindx,varname,timestrings(itime,:))};

   %% Determine the hyperslab indexing
   
   myinfo.diagn_file = fname;
   if (isfinite(copyindex)), myinfo.copyindex  = copyindex; end
   myinfo.levelindex = levelindx;
   myinfo.timeindex  = timeindx;
   [start, count]    = GetNCindices(myinfo, 'diagn', varname);
   
   % Extract field
   
   datmat = double(nc_varget(fname, varname, start, count));
   
   clf; 
   h = pcolor(xy.lonmat, xy.latmat, datmat);
   set(h,'LineStyle','none');
   
   title(plot_title,'Interpreter','none')
   h2 = colorbar('vert');
   set(get(h2,'YLabel'),'String',varinfo.units)

end



function varexist(filename, varnames)
%% We already know the file exists by this point.
% Lets check to make sure that file contains all needed variables.

nvars  = length(varnames);
gotone = ones(1,nvars);

for i = 1:nvars
   gotone(i) = nc_isvar(filename,varnames{i});
   if ( ~ gotone(i) )
      warning('\n%s is not a variable in %s\n',varnames{i},filename)
   end
end

if ~ all(gotone) 
   error('missing required variable ... exiting')
end



function xy = GetWRFlatlon(fname, varname);
%% Each of the WRF variables has a 'coordinate' attribute signifying which
% of the 6 possible lat/lon variables is appropriate.

coordinates{1} = sscanf(nc_attget(fname,varname,'coordinates'),'%s %*s');
coordinates{2} = sscanf(nc_attget(fname,varname,'coordinates'),'%*s %s');

latcoord       = find(strncmp('XLAT', coordinates, length('XLAT')) > 0);
loncoord       = find(strncmp('XLON', coordinates, length('XLON')) > 0);
xy.latmat      = nc_varget(fname, coordinates{latcoord});
xy.lonmat      = nc_varget(fname, coordinates{loncoord});
xy.latunits    = nc_attget(fname, coordinates{latcoord},'units');
xy.lonunits    = nc_attget(fname, coordinates{latcoord},'units');

inds = (xy.lonmat < 0); % Convert to 0,360 to minimize dateline probs.
xy.lonmat(inds) = xy.lonmat(inds) + 360.0;



function levels = GetWRFlevels(fname, varname);
%% Return the appropriate vertical indices.

varinfo = nc_getvarinfo(fname,varname);
lvlind  = find(strncmp('bottom_top',varinfo.Dimension,length('bottom_top')));
nlevels = varinfo.Size(lvlind);
levels  = 1:nlevels;



function times = GetWRFtimes(fname, varname);
%% Return the appropriate time coordinates.
%  DART files use 'time' as a real array
%  WRF  files use 'Time' as a character string

varinfo = nc_getvarinfo(fname,varname);
timeind = find(strncmp('time',varinfo.Dimension,length('time')));

if (isempty(timeind)) % WRF flavor Time variable.
    timeind = find(strncmp('Time',varinfo.Dimension,length('Time')));
    if (isempty(timeind)), error('%s has no time information. Dying.'); end
    timestr  = nc_varget(fname, 'Times');  % Note the plural ... seeth.
    timebase = sscanf(timestr,'%d%*c%d%*c%d%*c%d%*c%d%*c%d'); % YYYY MM DD HH MM SS
    times    = datenum(timebase');

else
    times       = nc_varget(fname, 'time');
    calendar    = nc_attget(fname, 'time','calendar');
    timeunits   = nc_attget(fname, 'time','units');
    timebase    = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
    timeorigin  = datenum(timebase(1),timebase(2),timebase(3));
    times       = times + timeorigin;
 
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
