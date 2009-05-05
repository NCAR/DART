function obsstruct = plot_obs_netcdf(fname, ObsTypeString, region, ObsString, ...
                                     QCString, maxQC, verbose)
%
% fname         = 'obs_sequence_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';
% region        = [0 360 -90 90 -Inf Inf];
% ObsString     = 'NCEP BUFR observation';
% QCString      = 'DART quality control';
% maxQC         = 2;
% verbose       = 1;   % anything > 0 == 'true'
%
% bob = plot_obs_netcdf(fname, ObsTypeString, region, ObsString, QCString, maxQC, verbose);

% record the user input

obsstruct.fname         = fname;
obsstruct.ObsTypeString = ObsTypeString;
obsstruct.region        = region;
obsstruct.ObsString     = ObsString;
obsstruct.QCString      = QCString;
obsstruct.maxQC         = maxQC;
obsstruct.verbose       = verbose;

% get going

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

my_types   = unique(obs_type);  % only ones in the file, actually.
timeunits  = nc_attget(fname,'time','units');
timerange  = nc_attget(fname,'time','valid_range');
calendar   = nc_attget(fname,'time','calendar');
timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
timestring = datestr(timerange + timeorigin);

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

% Find observations of the correct type.

mytypeind = get_copy_index(fname, ObsString);
myind     = strmatch(ObsTypeString,ObsTypeStrings);
inds      = find(obs_type == myind);
mylocs    = loc(inds,:);
myobs     = obs(inds,mytypeind);

if ~ isempty(QCString)
   myQCind = get_qc_index(fname,  QCString);
   myqc    = qc(inds,myQCind);
else
   myqc    = [];
end

% geographic subset if needed

inds = locations_in_region(mylocs,region);

obsstruct.lons = mylocs(inds,1);
obsstruct.lats = mylocs(inds,2);
obsstruct.z    = mylocs(inds,3);
obsstruct.obs  =  myobs(inds);

if ( ~ isempty(myqc))
obsstruct.qc   =   myqc(inds);
end

% It might be great to have a histogram of the observations with particular QC
% values

% subset based on qc value

if ( (~ isempty(myqc)) & ( ~ isempty(maxQC)) )

   inds = find(obsstruct.qc > maxQC);
   disp(sprintf('Removing %d obs with a %s value greater than %f', ...
                length(inds),QCString,maxQC))

   inds = find(obsstruct.qc <= maxQC);

   bob = obsstruct.lons(inds); obsstruct.lons = bob;
   bob = obsstruct.lats(inds); obsstruct.lats = bob;
   bob = obsstruct.z(   inds); obsstruct.z    = bob;
   bob = obsstruct.obs( inds); obsstruct.obs  = bob;

end

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
% axis([xmin xmax -Inf Inf])
worldmap;
colorbar;

str1 = sprintf('%s level (%.2f - %.2f)',ObsTypeString,zmin,zmax);
str2 = sprintf('(%d locations)',length(obsstruct.obs));
str3 = sprintf('%s - %s',timestring(1,:),timestring(2,:));

title( {str1, str3}, 'Interpreter','none','FontSize',18);
xlabel(str2,'FontSize',16)
