function obsstruct = read_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
   QCString, verbose)
%% read_obs_netcdf reads in the netcdf flavor observation sequence file
%                  and returns a subsetted structure.
%
% fname         = 'obs_epoch_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';   % or 'ALL' ...
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'NCEP BUFR observation';         % or 'ALL' ...
% QCString      = 'DART quality control';
% verbose       = 1;   % anything > 0 == 'true'
%
% obs = read_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, verbose);
%
% The return variable 'obs' is a structure. As an example ...
%
%           fname: 'obs_epoch_001.nc'
%   ObsTypeString: 'RADIOSONDE_U_WIND_COMPONENT'
%          region: [0 360 -90 90 -Inf Inf]
%      CopyString: 'NCEP BUFR observation'
%        QCString: 'DART quality control'
%         verbose: 1
%      timestring: [2x20 char]
%            lons: [2343x1 double]
%            lats: [2343x1 double]
%               z: [2343x1 double]
%             obs: [2343x1 double]
%            Ztyp: [2343x1 int32]
%            keys: [2343x1 int32]
%            time: [2343x1 double]
%              qc: [2343x1 int32]
%
%--------------------------------------------------
% EXAMPLE 1: plotting the locations without regard to observation value.
%--------------------------------------------------
% fname         = 'obs_epoch_001.nc';
% ObsTypeString = 'ALL';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'observation';
% QCString      = 'DART quality control ';
% verbose       = 1;   % anything > 0 == 'true'
%
% bob = read_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, verbose);
% plot(bob.lons,bob.lats,'*')
% continents('light');

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

%% record the user input

obsstruct.fname         = fname;
obsstruct.ObsTypeString = ObsTypeString;
obsstruct.region        = region;

%%
switch lower(CopyString)
   case 'all'
      obsstruct.CopyString = cellstr(ncread(fname,'CopyMetaData')');
   otherwise
      obsstruct.CopyString = CopyString;
end
obsstruct.QCString      = QCString;
obsstruct.verbose       = verbose;

%% get going

ObsTypes       = ncread(fname,'ObsTypes');
ObsTypeStrings = cellstr(ncread(fname,'ObsTypesMetaData')');
CopyStrings    = cellstr(ncread(fname,'CopyMetaData')');
QCStrings      = cellstr(ncread(fname,'QCMetaData')');

t              = ncread(fname,'time');
obs_type       = ncread(fname,'obs_type');
obs_keys       = ncread(fname,'obs_keys');

% FIXME ... if which_vert exists, use it.
% cartesian models do not have it
% if it doesn't exist - they all have the same z_type ?
if (nc_var_exists(fname,'which_vert'))
    z_type = ncread(fname,'which_vert');
else
    z_type = ones(size(obs_keys));
end

% ncread does not recognize the 'missing_value' attribute,
% it recognizes '_FillValue' ... sheesh

loc            = apply_missing(fname,'location');
obs            = apply_missing(fname,'observations');
qc             = apply_missing(fname,'qc');

my_types       = unique(obs_type);  % only ones in the file, actually.
timeunits      = nc_read_att(fname,'time','units');
timerange      = nc_read_att(fname,'time','valid_range');
calendar       = nc_read_att(fname,'time','calendar');
timebase       = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin     = datenum(timebase(1),timebase(2),timebase(3));
timestring     = datestr(timerange + timeorigin);
t              = t + timeorigin;

obsstruct.timestring = timestring;

%% Echo summary if requested

if ( verbose > 0 )
   for i = 1:length(my_types)
      obtype = my_types(i);
      inds   = find(obs_type == obtype);
      myz    = loc(inds,3);

      fprintf('N = %6d %32s (type %3d) tween levels %.2f and %.2f\n', ...
         length(inds), ObsTypeStrings{obtype}, obtype, ...
         unique(min(myz)), unique(max(myz)))
   end
end

%% Find copies of the correct type.
%  If 'ALL' is requested ... do not subset.

ncopies = nc_dim_info(fname,'copy');
switch lower(CopyString)
   case 'all'
      mytypeind = 1:ncopies;
   otherwise
      mytypeind = get_copy_index(fname, CopyString, 'context', 'CopyString');
end

%% Find observations of the correct type.
%  If 'ALL' is requested ... do not subset.

switch lower(ObsTypeString)
   case 'all'
      inds      = 1:size(obs,1);

   otherwise % subset the ONE desired observation type
      myind = strcmp(ObsTypeString, ObsTypeStrings);

      if ( any(myind) )
         myind = find(myind > 0,1); % find first instance of ...
         inds  = find(obs_type == myind);
      else
         fprintf('FYI - no %s observations ...\n',obsstruct.ObsTypeString)
         inds = [];
      end
end

myobs  =      obs(inds,mytypeind);
mylocs =      loc(inds,:);
mykeys = obs_keys(inds);
myztyp =   z_type(inds);
mytime =        t(inds);

%% Find desired QC values of those observations

if ~ isempty(QCString)
   myQCind = get_qc_index(fname,  QCString, 'QCString');
   myqc    = qc(inds,myQCind);
else
   myqc    = [];
end

%% geographic subset if needed

inds = locations_in_region(mylocs,region);

obsstruct.numobs = length(inds);
obsstruct.lons = mylocs(inds,1);
obsstruct.lats = mylocs(inds,2);
obsstruct.z    = mylocs(inds,3);
obsstruct.obs  =  myobs(inds,:);
obsstruct.Ztyp = myztyp(inds);
obsstruct.keys = mykeys(inds);
obsstruct.time = mytime(inds);
obsstruct.qc   = [];

if ~ isempty(myqc)
   obsstruct.qc = myqc(inds);
end

%% Try to determine if large numbers of the Z coordinate are up or down ...
%  The observation sequences don't actually have the units as part of the
%  metadata. That would be 'nice'.

ztypes = unique(obsstruct.Ztyp);

obsstruct.numZtypes = length(ztypes);

for itype = 1:obsstruct.numZtypes

   if (     ztypes(itype) == -2 )   % VERTISUNDEF
      obsstruct.Zpositivedir = 'up';
      obsstruct.Zunits       = 'unknown';

   elseif ( ztypes(itype) == -1 )   % VERTISSURFACE
      obsstruct.Zpositivedir = 'up';
      obsstruct.Zunits       = 'surface';

   elseif ( ztypes(itype) ==  1 )   % VERTISLEVEL
      % for example, cam level 1 is at the top of the atmosphere
      obsstruct.Zpositivedir = 'down';
      obsstruct.Zunits       = 'model level';

   elseif ( ztypes(itype) ==  2 )   % VERTISPRESSURE
      obsstruct.Zpositivedir = 'down';
      obsstruct.Zunits       = 'pressure';

   elseif ( ztypes(itype) ==  3 )   % VERTISHEIGHT
      %% here is the troublemaker.
      obsstruct.Zpositivedir = 'up';
      obsstruct.Zunits       = 'height';

      % If the observations are from the World Ocean Database, they
      % are probably depths. large positive numbers are near the
      % center of the earth. The thing is, the only reliable way to know
      % is to check the FIRST copy string.

      inds = strncmp(CopyStrings{1},'WOD obs',7);
      if ( inds > 0 )
         obsstruct.Zpositivedir = 'down';
         obsstruct.Zunits       = 'depth';
      end

   else
      error('not a known vertical coordinate system')
   end

end


function hyperslab = apply_missing(fname,varname)
% ncread does not recognize the 'missing_value' attribute,
% it recognizes '_FillValue' ... sheesh

hyperslab = ncread(fname,varname)';
missing   = nc_read_att(fname,varname,'missing_value');

if ~isempty(missing)
    hyperslab(hyperslab==missing) = NaN;
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
