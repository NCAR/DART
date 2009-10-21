function obsstruct = read_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
                                     QCString, maxQC, verbose)
%% read_obs_netcdf reads in the netcdf flavor observation sequence file
%                  and returns a subsetted structure.
%
% fname         = 'obs_sequence_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';   % or 'ALL' ...
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'NCEP BUFR observation';
% QCString      = 'DART quality control';
% maxQC         = 2;
% verbose       = 1;   % anything > 0 == 'true'
%
% obs = read_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, maxQC, verbose);
%
% The return variable 'obs' is a structure. As an example ...
%
%           fname: 'obs_sequence_001.nc'
%   ObsTypeString: 'RADIOSONDE_U_WIND_COMPONENT'
%          region: [0 360 -90 90 -Inf Inf]
%      CopyString: 'NCEP BUFR observation'
%        QCString: 'DART quality control'
%           maxQC: 2
%         verbose: 1
%      timestring: [2x20 char]
%            lons: [2343x1 double]
%            lats: [2343x1 double]
%               z: [2343x1 double]
%             obs: [2343x1 double]
%            Ztyp: [2343x1 double]
%        numbadqc: 993
%              qc: [2343x1 double]
%          badobs: [1x1 struct]

%% Data Assimilation Research Testbed -- DART
%  Copyright 2004-2007, Data Assimilation Research Section
%  University Corporation for Atmospheric Research
%  Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
%  <next few lines under version control, do not edit>
%  $URL$
%  $Id$
%  $Revision$
%  $Date$

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

%% record the user input

obsstruct.fname         = fname;
obsstruct.ObsTypeString = ObsTypeString;
obsstruct.region        = region;
obsstruct.CopyString    = CopyString;
obsstruct.QCString      = QCString;
obsstruct.maxQC         = maxQC;
obsstruct.verbose       = verbose;

%% get going

ObsTypes       = nc_varget(fname,'ObsTypes');
ObsTypeStrings = nc_varget(fname,'ObsTypesMetaData');
CopyStrings    = nc_varget(fname,'CopyMetaData');
QCStrings      = nc_varget(fname,'QCMetaData');

t              = nc_varget(fname,'time');
obs_type       = nc_varget(fname,'obs_type');
z_type         = nc_varget(fname,'which_vert');

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

obsstruct.timestring = timestring;

%% Echo summary if requested

if ( verbose > 0 ) 
   for i = 1:length(my_types)
      obtype = my_types(i);
      inds   = find(obs_type == obtype);
      myz    = loc(inds,3);
     
      fprintf('N = %6d %s obs (type %3d) between levels %.2f and %.2f\n', ...
               length(inds), ObsTypeStrings(obtype,:), obtype, ...
               unique(min(myz)), unique(max(myz)))
   end

%  uniquelevels = unique(loc(:,3));

%  for i = 1:length(uniquelevels)
%     mylevel = uniquelevels(i);
%     inds    = find(loc(:,3) == mylevel);
%     disp(sprintf('level %2d %f has %d observations',i,mylevel,length(inds)))
%  end

end

%% Find observations of the correct type.
%  If 'ALL' is requested ... do not subset.

mytypeind = get_copy_index(fname, CopyString);

switch lower(ObsTypeString)
   case 'all'
      inds      = 1:size(obs,1);

   otherwise % subset the desired observation type
      myind     = strmatch(ObsTypeString, ObsTypeStrings);
      if ( isempty(myind) ) 
         error('no %s observations ... stopping',obsstruct.ObsTypeString) 
      end
      inds      = find(obs_type == myind);
end

mylocs    = loc(inds,:);
myobs     = obs(inds,mytypeind);

%% Find desired QC values of those observations

if ~ isempty(QCString)
   myQCind = get_qc_index(fname,  QCString);
   myqc    = qc(inds,myQCind);
else
   myqc    = [];
end

%% geographic subset if needed

inds = locations_in_region(mylocs,region);

obsstruct.lons = mylocs(inds,1);
obsstruct.lats = mylocs(inds,2);
obsstruct.z    = mylocs(inds,3);
obsstruct.obs  =  myobs(inds);
obsstruct.Ztyp = z_type(inds);
obsstruct.qc   = [];
obsstruct.numbadqc = 0;

if ~ isempty(myqc)
   obsstruct.qc = myqc(inds);
end

%% subset based on qc value
%  

if ( (~ isempty(myqc)) && (~ isempty(maxQC)) )

   inds = find(obsstruct.qc > maxQC);

   obsstruct.numbadqc = length(inds);
   
   if (~isempty(inds))
       badobs.lons = obsstruct.lons(inds);
       badobs.lats = obsstruct.lats(inds);
       badobs.Ztyp = obsstruct.Ztyp(inds);
       badobs.z    = obsstruct.z(   inds);
       badobs.obs  = obsstruct.obs(inds);
       badobs.qc   = obsstruct.qc(inds);       
   end
   
   fprintf('Removing %d obs with a %s value greater than %f\n', ...
                length(inds),QCString,maxQC)

   inds = find(obsstruct.qc <= maxQC);

   bob = obsstruct.lons(inds); obsstruct.lons = bob;
   bob = obsstruct.lats(inds); obsstruct.lats = bob;
   bob = obsstruct.Ztyp(inds); obsstruct.Ztyp = bob;
   bob = obsstruct.z(   inds); obsstruct.z    = bob;
   bob = obsstruct.obs( inds); obsstruct.obs  = bob;
   bob = obsstruct.qc(  inds); obsstruct.qc   = bob;

end

if ( exist('badobs','var') )
   obsstruct.badobs = badobs;
end

