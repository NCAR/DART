function h = check_feom_profile(fname, vertex_list, varargin) 

% Simple Example:
%
% ncks -A -v longitudes,latitudes,depths Prior_Diag.nc Difference.nc
%
% fname = 'Difference.nc';
% vertex_list = 'temp_cell_7630_vertices.txt';
% h = check_feom_profile(fname,vertex_list);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

default_varname    = 'salt';
default_range      = [];
default_timestep   = 1;
default_level      = -1;
default_markersize =  8;
default_copystring = 'ensemble mean';
default_verbosity  = 0;  % 0 is silent
p = inputParser;

addRequired(  p, 'fname'       , @ischar);
addRequired(  p, 'vertex_list' , @ischar);
addParameter(p, 'varname'     , default_varname   , @ischar);
addParameter(p, 'range'       , default_range     , @isnumeric);
addParameter(p, 'time'        , default_timestep  , @isnumeric);
addParameter(p, 'level'       , default_level     , @isnumeric);
addParameter(p, 'markersize'  , default_markersize, @isnumeric);
addParameter(p, 'copystring'  , default_copystring, @ischar);
addParameter(p, 'verbosity'   , default_verbosity , @isnumeric);
parse(p, fname, vertex_list, varargin{:});

% if you want to echo the input
if (p.Results.verbosity ~= 0)
   fprintf('fname      : %s\n',     p.Results.fname)
   fprintf('vertex_list: %s\n',     p.Results.vertex_list)
   fprintf('varname    : %s\n',     p.Results.varname)
   fprintf('range      : %f %f \n', p.Results.range)
   fprintf('timestep   : %d\n',     p.Results.time)
   fprintf('level      : %d\n',     p.Results.level)
   fprintf('markersize : %d\n',     p.Results.markersize)
   fprintf('copystring : %s\n',     p.Results.copystring)
end

% do some error-checking on the input

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

if (exist(vertex_list,'file') ~= 2)
   error('%s does not exist.',vertex_list)
end

% gather the data for plotting

lons       = ncread(fname,'longitudes');
lats       = ncread(fname,'latitudes');
depths     = ncread(fname,'depths');
nVertices  = numel(lons);

inds = load(vertex_list,'-ascii');
inds = mod(inds-1,nVertices)+1;

copyindex = get_copy_index(fname, p.Results.copystring);
slab = get_hyperslab('fname',fname, ...
                  'varname',p.Results.varname, ...
                  'copyindex',copyindex, ...
                  'tindex',p.Results.time);

x = lons(inds);
y = lats(inds);
z = depths(inds);
c = slab(inds);
h = semilogy(c,z,'-*');

set(gca,'FontSize',20,'YDir','reverse')

xlabel(p.Results.varname)
ylabel('depth (m)')
string1 = sprintf('lat/lon %f, %f',y(1),x(1));
title({string1,fname},'Interpreter','none')
grid('on')
fprintf('%s\n',string1)

