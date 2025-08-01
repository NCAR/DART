function append_seg_data

% Change the path 'filepath' to netcdf files. 
% PARAM_FILE : Original 'parameters_dis_seg.nc' file
% CSV_TABLE  : Includes segment data from the shapefile 'river_midpoints.csv'
% NTWRK_FILE : Similar to 'parameters_dis_seg.nc' but with appended variables 

% The following variables are created and appended
% to the parameters_dis_seg.nc geometry file:
%
% - seg_lats      : Latitude of the individual segments (at midpoint?)
% - seg_lons      : Longitude of the individual segments 
% - poi_gauges    : USGS ID for the available gauges  
% - num_up_links  : Number of upstream segments at each reach
% - fromIndsStart : Index of first upstream segment
% - fromIndsEnd   : Index of last upstream segment
% - fromIndices   : Index into link/segment variables (upstream indices)
% 
% DART Software - Copyright UCAR. Open source software, provided "as is" 
% Webpage: https://dart.ucar.edu/

filepath = '/Users/gharamti/Documents/MATLAB/Derecho/HydroDART/pywatershed/'; 
table_f  = 'river_midpoints.csv';
source_f = 'parameters_dis_seg.nc';
dest_f   = 'parameters_dis_seg_app.nc'; 

NTWRK_FILE = strcat(filepath, dest_f);     % Appended nc file
PARAM_FILE = strcat(filepath, source_f);   % Original params file
CSV_TABLE  = strcat(filepath, table_f);    % CSV table file

if ~isfile(PARAM_FILE) || ~isfile(CSV_TABLE)
    error([ 'Either the PARAM_FILE "' PARAM_FILE '" or the CSV_TABLE "' ...
          CSV_TABLE '" do not exist.' ])
end


%% Get lon from table (ArcGIS)
data = readtable(table_f);

inds = table2array(data(:, 5));
lats = table2array(data(:, 6));
lons = table2array(data(:, 7));

[~, isort] = sort(inds);

lats = lats(isort);
lons = lons(isort);

n_char = 15;

%% Compute Routelink style variables
ncPARAM = ncinfo(PARAM_FILE);

% read dimensions
n_seg = ncPARAM.Dimensions(1).Length;
gages = ncPARAM.Dimensions(2).Length;

% read variables
to_index   = double(ncread(PARAM_FILE, 'tosegment'));
poi_gauges = char(pad(ncread(PARAM_FILE, 'poi_gage_id'), n_char, 'left')); 

num_up_links  = zeros(1, n_seg);
upstream_seg  = struct;
fromIndsStart = zeros(1, n_seg);
fromIndsEnd   = zeros(1, n_seg);

k = 0;

for i = 1:n_seg
    
    % Upstream segments of i
    tmp = find(to_index == i)';

    % Find how many links are above me: 1, 2, 3, 4, ...
    % and get their indices. 
    num_up_links(i) = length(tmp);
    upstream_seg(i).links(:) = tmp;

    % If `fromIndsStart` and `fromIndsEnd` are the same 
    % and non-zero, that means there is only one link above me
    % num_up_links   = fromIndsEnd - fromIndsStart + 1; 
    % num_up_links(no_up_links) = 0;
    if ~isempty(tmp)
        fromIndsStart(i) = k+1;

        % Link i has upstream links
        for j = 1:num_up_links(i)
            k = k + 1;

            % Index into the state
            fromIndices(k) = tmp(j); %#ok
        end
        fromIndsEnd(i)   = k;
    end
end
n_fi = length(fromIndices);


%% Write NETCDF file
if isfile(NTWRK_FILE)
    warning([ 'A previosuly created "' dest_f '" file already exists. Deleting it ..' ])
    delete(NTWRK_FILE)
end

ncid = netcdf.create(NTWRK_FILE, 'NETCDF4'); % Full NetCDF-4 format

% Copy dimensions
for i = 1:length(ncPARAM.Dimensions)
    dim_name   = ncPARAM.Dimensions(i).Name;
    dim_length = ncPARAM.Dimensions(i).Length;
    
    % Use low-level netCDF function to define dimensions properly
    dim_ids.(dim_name) = netcdf.defDim(ncid, dim_name, dim_length);
end

% Define variables without data
var_ids = struct(); % Store variable IDs for reference
for i = 1:length(ncPARAM.Variables)
    var_name = ncPARAM.Variables(i).Name;
    var_dims = {ncPARAM.Variables(i).Dimensions.Name};  % Get dimension names
    var_type = ncPARAM.Variables(i).Datatype;

    % Convert dimension names to corresponding dimension IDs
    dim_id_list = [];
    for j = 1:length(var_dims)
        dim_id_list = [dim_id_list, dim_ids.(var_dims{j})]; %#ok
    end

    % Define the variable in the new NetCDF file
    var_ids.(var_name) = netcdf.defVar(ncid, var_name, var_type, dim_id_list);
end

% Exit definition mode
netcdf.endDef(ncid);

% Copy variable data
for i = 1:length(ncPARAM.Variables)
    var_name = ncPARAM.Variables(i).Name;
    data     = ncread(PARAM_FILE, var_name);
    
    % Write data to the new file
    netcdf.putVar(ncid, var_ids.(var_name), data);
end

nccreate(NTWRK_FILE, 'fromIndsStart', 'Dimensions', {'nsegment' , n_seg}, 'DataType', 'int64')
nccreate(NTWRK_FILE, 'fromIndsEnd'  , 'Dimensions', {'nsegment' , n_seg}, 'DataType', 'int64')
nccreate(NTWRK_FILE, 'fromIndices'  , 'Dimensions', {'index'    , n_fi }, 'DataType', 'int64')
nccreate(NTWRK_FILE, 'num_up_links' , 'Dimensions', {'nsegment' , n_seg}, 'DataType', 'int64')
nccreate(NTWRK_FILE, 'seg_lats'     , 'Dimensions', {'nsegment' , n_seg}, 'DataType', 'double')
nccreate(NTWRK_FILE, 'seg_lons'     , 'Dimensions', {'nsegment' , n_seg}, 'DataType', 'double')
nccreate(NTWRK_FILE, 'poi_gauges'   , 'Dimensions', {'IDLength' , n_char, 'npoigages', gages}, 'DataType', 'char')

ncwrite(NTWRK_FILE, 'fromIndsStart', fromIndsStart)
ncwrite(NTWRK_FILE, 'fromIndsEnd'  , fromIndsEnd)
ncwrite(NTWRK_FILE, 'fromIndices'  , fromIndices)
ncwrite(NTWRK_FILE, 'num_up_links' , num_up_links)
ncwrite(NTWRK_FILE, 'seg_lats'     , lats)
ncwrite(NTWRK_FILE, 'seg_lons'     , lons)
ncwrite(NTWRK_FILE, 'poi_gauges'   , poi_gauges')

netcdf.close(ncid); 
