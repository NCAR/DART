function plotdat = read_obsdiag_staticdata(fname, copy)
%% read the static data from the netCDF output of threed_sphere/obs_diag

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

plotdat = struct('fname',fname,'copystring',copy);

plotdat.bincenters    = ncread(fname,'time');
plotdat.binedges      = ncread(fname,'time_bounds');
plotdat.mlevel        = local_ncread(fname,'mlevel');
plotdat.plevel        = local_ncread(fname,'plevel');
plotdat.plevel_edges  = local_ncread(fname,'plevel_edges');
plotdat.hlevel        = local_ncread(fname,'hlevel');
plotdat.hlevel_edges  = local_ncread(fname,'hlevel_edges');
[plotdat.ncopies, ~]  = nc_dim_info(fname,'copy');
[plotdat.nregions, ~] = nc_dim_info(fname,'region');
plotdat.region_names  = strtrim(ncread(fname,'region_names')');

plotdat.dimensionality = nc_read_att(fname, '/', 'LocationRank');
plotdat.binseparation = nc_read_att(fname, '/', 'bin_separation');
plotdat.binwidth      = nc_read_att(fname, '/', 'bin_width');
plotdat.lonlim1       = nc_read_att(fname, '/', 'lonlim1');
plotdat.lonlim2       = nc_read_att(fname, '/', 'lonlim2');
plotdat.latlim1       = nc_read_att(fname, '/', 'latlim1');
plotdat.latlim2       = nc_read_att(fname, '/', 'latlim2');
plotdat.biasconv      = nc_read_att(fname, '/', 'bias_convention');

plotdat.copyindex     = get_copy_index(fname,copy);
plotdat.rmseindex     = get_copy_index(fname,'rmse');
plotdat.Npossindex    = get_copy_index(fname,'Nposs');
plotdat.Nusedindex    = get_copy_index(fname,'Nused');
plotdat.NQC0index     = get_copy_index(fname,'N_DARTqc_0');
plotdat.NQC1index     = get_copy_index(fname,'N_DARTqc_1');
plotdat.NQC2index     = get_copy_index(fname,'N_DARTqc_2');
plotdat.NQC3index     = get_copy_index(fname,'N_DARTqc_3');
plotdat.NQC4index     = get_copy_index(fname,'N_DARTqc_4');
plotdat.NQC5index     = get_copy_index(fname,'N_DARTqc_5');
plotdat.NQC6index     = get_copy_index(fname,'N_DARTqc_6');
plotdat.NQC7index     = get_copy_index(fname,'N_DARTqc_7');
plotdat.NQC8index     = get_copy_index(fname,'N_DARTqc_8','fatal',false);

% Coordinate between time types and dates

time_to_skip     = nc_read_att(fname, '/', 'time_to_skip');
timeunits        = nc_read_att(fname,'time','units');
timebase         = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin       = datenum(timebase(1),timebase(2),timebase(3));
if ( isempty(time_to_skip) == 1)
    iskip = 0;
elseif ( numel(time_to_skip) == 6)
    skip_seconds = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
    iskip        = time_to_skip(3) + skip_seconds/86400;
else
    error('time_to_skip variable has unusual length. Should be either 0 or 6.')
end

% Set up a structure to use for plotting

plotdat.bincenters = plotdat.bincenters + timeorigin;
plotdat.binedges   = plotdat.binedges   + timeorigin;
plotdat.Nbins      = length(plotdat.bincenters);
plotdat.toff       = plotdat.bincenters(1) + iskip;

%=====================================================================

function value = local_ncread(fname,varname)
%% If the variable exists in the file, return the contents of the variable.
% if the variable does not exist, return empty value instead of error-ing
% out.

[variable_present, ~] = nc_var_exists(fname,varname);
if (variable_present)
    value = ncread(fname,varname);
else
    value = [];
end
% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
