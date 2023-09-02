function even_sphere_compressed(nprofiles, varargin)

% Generate approximately evenly-distributed profiles on a sphere using Golden Section spiral algorithm
%
% This creates input to create_obs_sequence that will generate a set of radiosonde
% profiles of T,U,V, that are approximately evenly-spaced on a sphere.
% A subset of obs in a lon-lat box may be selected, but the spacing will be derived 
% from nprofiles *on the whole sphere*.
% An optional argument 'fill_obs' will trigger the creation of synthetic
% observations complete with (silly) observation values.
%
% The vertical coordinate system is 'pressure'.
% The pressure levels themselves may be input.
% The observation error variances for each level must then also be input.
% There are separate error variances for T and [U,V].
%
% Example 1: 30 profiles, 14 default levels and obs error variances
% nprofiles   = 30;
% even_sphere_subdom(nprofiles)
%
% Example 2: 30 profiles at 6 specified levels, 
%            uses the default obs error variances for each layer
% nprofiles   = 30;
% levels      = [1000  850  500  300  200  100];
% even_sphere_subdom(nprofiles, 'levels', levels)
%
% Example 3: 30 profiles at 6 specified levels, 
%            specify the error variances for each layer 
% nprofiles   = 30;
% levels      = [1000  850  500  300  200  100];
% T_error_var = [1.44 0.64 0.64 0.81 1.44 0.64];
% W_error_var = [1.96 2.25 4.41 9.00 7.29 4.41];
% even_sphere_subdom(nprofiles, 'levels', levels, ...
%             'T_error_var', T_error_var, 'W_error_var', W_error_var)
%
% Example 4: restrict profiles to a horizontal subdomain (which spans the prime meridian).
% nprofiles   = 2000;  only some of these will be written to the file.
% lon_left = 300.; left boundary
% lon_right =  30.; right boundary
% lat_south = -20.; southern boundary
% lat_north =  20.; northern boundary

% 

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------

p = inputParser;

addRequired(p,'nprofiles',@isnumeric);

% obs types to create at each location.
obs_types = {'RADIOSONDE_TEMPERATURE','RADIOSONDE_U_WIND_COMPONENT','RADIOSONDE_V_WIND_COMPONENT'};
num_types = size(obs_types,2);

% These are the *mandatory pressure levels* defined in the
% AMS glossary https://glossary.ametsoc.org/wiki/Mandatory_level
% The error variances come from obs_converters/obs_error/ncep_obs_err_mod.f90.

default_levels      = [1000  925  850  700  500  400  300   250  200  150  100   70   50   30   20   10    7    5    3    2    1];
default_T_error_var = [1.44 1.00 0.64 0.64 0.64 0.64 0.81  1.44 1.44 1.00 0.64 0.64 0.81 1.00 1.69 2.25 2.25 2.25 2.25 2.25 2.25];
default_W_error_var = [1.96 2.25 2.25 2.56 4.41 6.76 9.00 10.24 7.29 5.76 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41];
default_fill        = false; % For diagnostic test need a null data value and a null qc value
default_nlevels     = length(default_levels);
default_YMD         = '2017-12-25';
default_nlon         = 288;
default_nlat         = 192;
default_lon_left         =   0.;
default_lon_right         = 360.;
default_lat_south         = -90.;
default_lat_north         =  90.;

if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'T_error_var', default_T_error_var, @isnumeric);
    addParameter(p,'W_error_var', default_W_error_var, @isnumeric);
    addParameter(p,'levels',      default_levels,      @isnumeric);
    addParameter(p,'nlevels',     default_nlevels,     @isnumeric);
    addParameter(p,'fill_obs',    default_fill,        @islogical);
    addParameter(p,'YMD',         default_YMD,         @ischar);
    addParameter(p,'nlon',        default_nlon,        @isnumeric);
    addParameter(p,'nlat',        default_nlat,        @isnumeric);
    addParameter(p,'lon_left',        default_lon_left,        @isnumeric);
    addParameter(p,'lon_right',        default_lon_right,        @isnumeric);
    addParameter(p,'lat_south',        default_lat_south,        @isnumeric);
    addParameter(p,'lat_north',        default_lat_north,        @isnumeric);
    fprintf('called addParameters ')
else
    addParamValue(p,'T_error_var', default_T_error_var, @isnumeric); %#ok<NVREPL>
    addParamValue(p,'W_error_var', default_W_error_var, @isnumeric); %#ok<NVREPL>
    addParamValue(p,'levels',      default_levels,      @isnumeric); %#ok<NVREPL>
    addParamValue(p,'nlevels',     default_nlevels,     @isnumeric); %#ok<NVREPL>
    addParamValue(p,'fill_obs',    default_fill,        @islogical); %#ok<NVREPL>
    addParamValue(p,'YMD',         default_YMD,         @ischar);    %#ok<NVREPL>
    addParamValue(p,'nlon',        default_nlon,        @isnumeric); %#ok<NVREPL>
    addParamValue(p,'nlat',        default_nlat,        @isnumeric); %#ok<NVREPL>
    addParamValue(p,'lon_left',        default_lon_left,        @isnumeric);
    addParamValue(p,'lon_right',        default_lon_right,        @isnumeric);
    addParamValue(p,'lat_south',        default_lat_south,        @isnumeric);
    addParamValue(p,'lat_north',        default_lat_north,        @isnumeric);
    fprintf('called addParamValues ')
end

p.parse(nprofiles, varargin{:});

nlevels     = p.Results.nlevels;
levels      = p.Results.levels;
T_error_var = p.Results.T_error_var;
W_error_var = p.Results.W_error_var;
fill_obs    = p.Results.fill_obs;
yyyymmdd    = p.Results.YMD;
nlon        = p.Results.nlon;
nlat        = p.Results.nlat;
lon_left        = p.Results.lon_left;
lon_right        = p.Results.lon_right;
lat_south        = p.Results.lat_south;
lat_north        = p.Results.lat_north;

% FIXME ... check that all arrays of length 'level' are consistent.

% If the user specifies nlevels, use the first nlevels of the default.
% If the user also specifies 'levels', use the first nlevels of what they
% specify.

nlevels    = min(nlevels,length(levels));
levels     =       levels(1:nlevels);
T_error_var = T_error_var(1:nlevels);
W_error_var = W_error_var(1:nlevels);

levelstrings = sprintf('%.2f ',levels);
fprintf('vertical levels %s\n',levelstrings)

% Generate obs_sequence input for this problem

% If create_fixed_network_sequence is run, this date information will be overwritten

[year, month, day] = ymd(datetime(yyyymmdd));
obsdate   = datenum(yyyymmdd);
gregorian = datenum('1601-01-01');
dart_date = obsdate-gregorian;

fprintf('Using date/time as %d-%02d-%02d 00:00:00 Z\n',year, month, day)
fprintf('(days since 1601-01-01) : %f\n',dart_date)

fprintf('Creating %d profiles with %d levels for %d variables = %d observations ...\n', ...
         nprofiles,nlevels,num_types,nprofiles*nlevels*num_types)

% preallocate space for efficiency
x(       1:nprofiles) = 0;
y(       1:nprofiles) = 0;
z(       1:nprofiles) = 0;
lon(     1:nprofiles) = 0;   % radians
lat(     1:nprofiles) = 0;   % radians
lon_deg( 1:nprofiles) = 0;   % degrees
lat_deg( 1:nprofiles) = 0;   % degrees
plot_lon(1:nprofiles) = 0;   % radians, Only num_plot will be filled
plot_lat(1:nprofiles) = 0;   % radians, Only num_plot will be filled

% Calculate the cartesian and spherical locations.
% For the geometric and visually minded: see the README.rst.
inc = pi * (3 - sqrt(5));
off = 2 / nprofiles;
for k = 1:nprofiles
    y(k) = (k-1) * off - 1 + (off / 2);
    r = sqrt(1 - y(k)*y(k));
    phi = (k-1) * inc;
    x(k) = cos(phi) * r;
    z(k) = sin(phi) * r;
    lon(k) = atan2(y(k), x(k)) + pi;
    lat(k) = asin(z(k));
end

% Now convert latitude and longitude to degrees and restrict to the subdomain.

num_plot = 0;
for k = 1:nprofiles
    % Input is in degrees of latitude and longitude; convert from radians
    lon_deg(k) = rad2deg(lon(k));
    lat_deg(k) = rad2deg(lat(k));

    % Restrict to the subdomain.
    % ? Should the poles be excluded?
    if ( lat_deg(k) >= lat_south & lat_deg(k) <= lat_north) 
       % Does the subdomain span the prime meridian?
       if (lon_left < lon_right) 
          if (lon_deg(k) >= lon_left & lon_deg(k) <= lon_right) 
             num_plot = num_plot + 1 ;
             plot_lon(num_plot) = lon(k);
             plot_lat(num_plot) = lat(k);
          end
       else
          if (lon_deg(k) >= lon_left || lon_deg(k) <= lon_right)
             num_plot = num_plot + 1 ;
             plot_lon(num_plot) = lon(k);
             plot_lat(num_plot) = lat(k);
          end
       end
    %    fprintf('%d is in the box; %f , %f \n',k, lon_deg(k), lat_deg(k))
    end
end
% Total number of observations at single time is nlevels*nprofiles*[T,U,V,...],
num_obs = num_plot * nlevels * num_types

% text file for driving create_obs_sequence
outputfile = 'even_create_input';
fprintf('Creating output file "%s" ... ',outputfile)
fid = fopen(outputfile, 'w');

% Output the total number of observations
fprintf(fid, '%6i\n', num_obs);

if(fill_obs)
    fprintf(fid, '%2i\n', 1);
    fprintf(fid, '%2i\n', 1);
    fprintf(fid, '"observations"\n'); % Metadata for data copy
    fprintf(fid, '"Quality Control"\n'); % Metadata for qc
else
    fprintf(fid, '%2i\n', 0); % 0 copies of data
    fprintf(fid, '%2i\n', 0); % 0 QC fields
end


% Loop to create each observation
for hloc = 1:num_plot
    for vloc = 1:length(levels)
        for field = obs_types
            
            % 0 indicates that there is another observation;
            fprintf(fid, '%2i\n', 0);
            
            % Specify obs kind by string
            fprintf(fid, '%s\n',string(field));
            
            % Specify pressure as the vertical coordinate system
            fprintf(fid, '%2i\n', 2);
            
            % The vertical pressure value
            fprintf(fid, '%5i\n', levels(vloc));
            
            % longitude and latitude in degrees
            fprintf(fid, '%6.2f\n', lon_deg(hloc));
            fprintf(fid, '%6.2f\n', lat_deg(hloc));
            
            % Now the date and time
            fprintf(fid, '%5i %3i %3i %3i %2i %2i \n', year, month, day, 0, 0, 0);
            
            % Finally, the error variance
            switch string(field)
                case 'RADIOSONDE_TEMPERATURE'
                    fprintf(fid, '%6.2f\n', T_error_var(vloc));
                otherwise
                    fprintf(fid, '%6.2f\n', W_error_var(vloc));
            end
            
            % Need observation value and qc for testing
            if(fill_obs)
                fprintf(fid, '%2i\n', 1);  % observation value
                fprintf(fid, '%2i\n', 0);  % qc value
            end
        end
    end
end

% File name that results from running create_obs_sequence.
fprintf(fid, 'set_def.out\n');

% Done with output, close file
fclose(fid);
fprintf('done.\n')


%-------------------------------------------------------------------------------
%% Some plotting to visualize this observation set in context
%  Plot a regular grid for comparison.

% Model grid for comparison on the obs distribution plot
% T85          : nlon=256, nlat=128 points on A-grid
% FV ~1 degree : nlon=288, nlat=192
% 1 degree every 4th point in latitude and longitude would be (288*192/16=3456)

% plot(lon, lat, 'x');
plot(plot_lon(1:num_plot), plot_lat(1:num_plot), 'x');
title(sprintf('%d profile locations with %d levels (1 level shown)',num_plot,nlevels))
set(gca,'FontSize',20)
axis([0. 2*pi -pi/2 pi/2])

yticks = -pi/2 :   pi/nlat : pi/2;
xticks =     0 : 2*pi/nlon : 2*pi;

% set(gca,'YTick',yticks, 'YTickLabel',[], ...
%     'XTick',    xticks, 'XTickLabel',[])
set(gca, 'YTick', yticks, 'YTickLabel',[], ...
         'XTick', xticks, 'XTickLabel',[])
grid on

xlabel(sprintf('%d evenly-spaced longitudes',nlon))
ylabel(sprintf('%d evenly-spaced latitudes' ,nlat))

% axis tight
