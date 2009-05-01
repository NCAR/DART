function obsstruct = plot_obs_netcdf(fname, ObsTypeString, region, ObsString, QCString, verbose)
%
% fname         = 'obs_sequence_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';
% region        = [0 360 -90 90 -Inf Inf];
% ObsString     = 'NCEP BUFR observation';
% QCString      = 'DART quality control';
% verbose       = 1;   % anything > 0 == 'true'
%
% bob = plot_obs_netcdf(fname, ObsTypeString, region, ObsString, QCString, verbose);

ObsTypes       = nc_varget(fname,'ObsTypes');
ObsTypeStrings = nc_varget(fname,'ObsTypesMetaData');
CopyStrings    = nc_varget(fname,'CopyMetaData');
QCStrings      = nc_varget(fname,'QCMetaData');

t              = nc_varget(fname,'time');
obs_type       = nc_varget(fname,'obs_type');
z              = nc_varget(fname,'which_vert');

loc            = nc_varget(fname,'location');
obs            = nc_varget(fname,'observations');
qc             = nc_varget(fname,'qc');

my_types = unique(obs_type);  % only ones in the file, actually.

% Echo summary if requested

if ( verbose > 0 ) 
   for i = 1:length(my_types)
      obtype = my_types(i);
      inds   = find(obs_type == obtype);
      myz    = loc(inds,3);
     
      disp(sprintf('N = %6d %s obs (type %3d) between levels %.2f and %.2f', ...
               length(inds), ObsTypeStrings(obtype,:), obtype, ...
               unique(min(myz)), unique(max(myz))))
   end

%  uniquelevels = unique(loc(:,3));

%  for i = 1:length(uniquelevels)
%     mylevel = uniquelevels(i);
%     inds    = find(loc(:,3) == mylevel);
%     disp(sprintf('level %2d %f has %d observations',i,mylevel,length(inds)))
%  end

end

mytypeind = get_copy_index(fname, ObsString);
myQCind   =   get_qc_index(fname,  QCString);

% Find observations of a the correct type.

myind  = strmatch(ObsTypeString,ObsTypeStrings);
inds   = find(obs_type == myind);
mylocs = loc(inds,:);
myqc   =  qc(inds,myQCind);
myobs  = obs(inds,mytypeind);

% geographic subset if needed

inds = locations_in_region(mylocs,region);

obsstruct.lons = mylocs(inds,1);
obsstruct.lats = mylocs(inds,2);
obsstruct.z    = mylocs(inds,3);
obsstruct.obs  =  myobs(inds);
obsstruct.qc   =   myqc(inds);

xmin = min(region(1:2));
xmax = max(region(1:2));
ymin = min(region(3:4));
ymax = max(region(3:4));
zmin = min(obsstruct.z);
zmax = max(obsstruct.z);

y1 = 36;
yN = 10*y1;
x1 = min(obsstruct.obs);
xN = max(obsstruct.obs);
slope = (yN-y1)/(xN-x1);
b = y1 - slope*x1;

scalarray = obsstruct.obs*slope + b;

h = scatter(obsstruct.lons,obsstruct.lats,scalarray,obsstruct.obs,'d','filled');
axis image
axis([xmin xmax ymin ymax])
worldmap;
colorbar;
title(sprintf('level (%.2f - %.2f) %s (%d locations)',zmin,zmax, ...
       ObsTypeString,length(inds)), 'Interpreter','none');

