function h = plot_feom_state(fname, varargin) 

% Simple Example:
%
% fname = '/glade/scratch/thoar/FILTER/Prior_Diag.nc';
% h = plot_feom_state(fname,'markersize',10);
%
% More complicated Example:
%
% h = plot_feom_state(fname,'varname','temp','timestep', 1, 'copystring','ensemble spread') 

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

addRequired(  p, 'fname'     , @ischar);
addParameter(p, 'varname'   , default_varname   , @ischar);
addParameter(p, 'range'     , default_range     , @isnumeric);
addParameter(p, 'time'      , default_timestep  , @isnumeric);
addParameter(p, 'level'     , default_level     , @isnumeric);
addParameter(p, 'markersize', default_markersize, @isnumeric);
addParameter(p, 'copystring', default_copystring, @ischar);
addParameter(p, 'verbosity' , default_verbosity , @isnumeric);
parse(p, fname, varargin{:});

% if you want to echo the input
if (p.Results.verbosity ~= 0)
   fprintf('fname      : %s\n',     p.Results.fname)
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

if ( ~ isempty(p.Results.range) && numel(p.Results.range) ~= 2)
   error('If specified, range must have exactly two elements.')
else
   default_range = sort(p.Results.range);
end

% gather the data for plotting

lons       = ncread(fname,'longitudes');
lats       = ncread(fname,'latitudes');
depths     = ncread(fname,'depths');
node_table = ncread(fname,'node_table'); % necessary to subset by level

copyindex = get_copy_index(fname, p.Results.copystring);

slab = get_hyperslab('fname',fname, ...
                  'varname',p.Results.varname, ...
                  'copyindex',copyindex, ...
                  'tindex',p.Results.time);

% subset the data to a plottable amount, 
% 3M points seems to be more than Matlab can handle.
% If the node_table returns a negative value, that is not a 'wet' vertex.
% only 'wet' vertices are in the FeoM state.

if (p.Results.level > 0) 
   vertices = node_table(p.Results.level,:);
   inds     = vertices(vertices > 0);
else
   error('cannot plot entire state ... must implement geographic subsetting.')
end

x = lons(inds);
y = lats(inds);
z = depths(inds);
s = zeros(size(x)) + p.Results.markersize;
c = slab(inds);
h = scatter(x,y,s,c,'filled');

set(gca,'FontSize',20)

if ( ~ isempty(default_range))
   set(gca,'Clim',default_range)
end

% h = scatter3(x,y,z,s,c,'filled');
% zlabel('depth (meters)')

xlabel('longitude (degrees East)')
ylabel('latitude (degrees North)')
string1 = sprintf('%s at level %d (%.1f m)',p.Results.varname,p.Results.level,z(1));
title({string1,fname},'Interpreter','none')
colorbar;
grid('on')

