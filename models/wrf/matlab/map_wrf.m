function field = map_wrf(fname, varname, copystring, levelindx, timeindx )
%% map_wrf creates a map of a WRF field - please check what is supposed to
%  happen with the staggering.
%
% (rather slow) Example: 
% fname      = 'Prior_Diag.nc';
% varname    = 'U_d01';
% copystring = 'ensemble mean';
% levelindx  = 10;
% timeindx   = 1;
%
% map_wrf(fname, varname, copystring, levelindx, timeindx );

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( nargin == 3 )
   levelindx = 1;
   timeindx  = 1;
elseif (nargin == 4 ) 
   timeindx  = 1;
elseif (nargin ~= 5 ) 
   error('Wrong number of arguments ... must have 3,4, or 5')
end

if (exist(fname,'file') ~= 2) 
   error('%s does not exist',fname)
end

varexist(fname,{varname,'time','XLON_d01','XLAT_d01','level_d01', ...
       'TRUELAT1','TRUELAT2', 'CEN_LAT', 'CEN_LON', 'MAP_PROJ'})

copyindex = get_copy_index(fname,copystring);

%% Read map metadata
times    = nc_varget(fname,     'time'); Ntimes = size(times, 1);
xlon     = nc_varget(fname, 'XLON_d01');  we    = size( xlon, 2);
xlat     = nc_varget(fname, 'XLAT_d01');  sn    = size( xlat, 1);
level    = nc_varget(fname,'level_d01');  bt    = size(level, 1);

stdlat1  = nc_varget(fname, 'TRUELAT1');
stdlat2  = nc_varget(fname, 'TRUELAT2');
cen_lat  = nc_varget(fname,  'CEN_LAT');
cen_lon  = nc_varget(fname,  'CEN_LON');
mp       = nc_varget(fname, 'MAP_PROJ')+1;  % convert [0,N-1] to [1,N]
map_proj = nc_attget(fname, 'MAP_PROJ','units');
inds     = find((map_proj == ',') | (map_proj == '='));
map_proj(inds) = ' ';
[a b] = strread(map_proj,'%s %s','delimiter',' ');
map_proj = lower(b);

minlat = min(xlat(:)); maxlat = max(xlat(:));
minlon = min(xlon(:)); maxlon = max(xlon(:));

%% Set up map projection
%  may have to translate their projections to matlabs keywords

fprintf('map projection is %s\n',map_proj{mp})

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

% iso = [0.5:1:5];
% iso = [0.01:0.01:0.1];
% iso = [50:50:300];
% iso = [0.5:0.5:5];
% iso = [100:100:600];
% iso = [0.0001:0.0001:0.001];
% iso = [0.00001:0.00001:0.0001];
% iso = [0.00001:0.00001:0.0001];

for itime = timeindx

   plot_title = sprintf('%s %s levelindex %d time %f',fname,varname,levelindx,times(timeindx));

   %% Determine the hyperslab indexing
   
   myinfo.diagn_file = fname;
   myinfo.copyindex  = copyindex;
   myinfo.levelindex = levelindx;
   myinfo.timeindex  = timeindx;
   [start, count]    = GetNCindices(myinfo, 'diagn', varname);
   
   % Extract field
   
   stag_field = nc_varget(fname, varname, start, count);
   field      = zeros(sn,we);

   switch lower(varname)
     case {'u','u10_d01','u_d01'}
        for iy = 1:sn
        for ix = 1:we
           field(iy,ix) = (stag_field(iy,ix) + stag_field(iy,ix+1))/2.0;
        end
        end
     case {'v','v10_d01','v_d01'}
        for iy = 1:sn
        for ix = 1:we
           field(iy,ix) = (stag_field(iy,ix) + stag_field(iy+1,ix))/2.0;
        end
        end
     otherwise
        field = stag_field;
   end
   
   % Plot field
   
   % subplot(m,m,pane);
   
   axesm(map_proj{mp}, 'Origin', [0 cen_lon 0], ...
         'MapParallels',[stdlat1 stdlat2],...
         'MapLatLimit',[minlat maxlat], ...
         'MapLonLimit',[minlon maxlon]);
   % framem;
   
%  states = shaperead('usastatelo', 'UseGeoCoords', true);
%  geoshow([states.Lat],[states.Lon],'Color','black');

   landareas = shaperead('landareas', 'UseGeoCoords', true);
   geoshow([landareas.Lat],[landareas.Lon],'Color','black');
   
   if min(field(:)) ~= max(field(:))

%     disp('contouring (slow) ...')
%     [C h]=contourm(xlat,xlon,field) ;
%     h = clabelm(C,h,'labelspacing',288);  set(h,'Fontsize',9);
   
      h1 = surfm(xlat,xlon,field);
   
      title(plot_title,'Interpreter','none')
      colorbar('vert')
   else
      error('field is all uniform value')
   end

end



function varexist(filename, varnames)
%% We already know the file exists by this point.
% Lets check to make sure that file contains all needed variables.

nvars  = length(varnames);
gotone = ones(1,nvars);

for i = 1:nvars
   gotone(i) = nc_isvar(filename,varnames{i});
   if ( ~ gotone(i) )
      fprintf('\n%s is not a variable in %s\n',varnames{i},filename)
   end
end

if ~ all(gotone) 
   error('missing required variable ... exiting')
end

