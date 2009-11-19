function link_obs(fname,ObsTypeString)
% link_obs generates the 'brushable' observation plots.
%
% EXAMPLE 1:
% fname         = 'obs_sequence_013.nc';
% ObsTypeString = 'RADIOSONDE_TEMPERATURE';
% link_obs(fname,ObsTypeString)
%
% EXAMPLE 3:
% link_obs('obs_sequence_003.nc','RADIOSONDE_TEMPERATURE')

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

region        = [0 360 -90 90 -Inf Inf];
CopyString    = 'observation';
CopyString    = 'NCEP BUFR observation';
QCString      = 'DART quality control';
verbose       = 1;

obs = read_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, verbose);

obs.lonindex  = 1;
obs.latindex  = 2;
obs.zindex    = 3;
obs.obsindex  = 4;
obs.qcindex   = 5;
obs.keyindex  = 6;
obs.timeindex = 7;
obs.indindex  = 8;

global obsmat
obsmat = zeros(length(obs.lons),5);
obsmat(:,obs.lonindex ) = obs.lons;
obsmat(:,obs.latindex ) = obs.lats;
obsmat(:,obs.zindex   ) = obs.z;
obsmat(:,obs.obsindex ) = obs.obs;
obsmat(:,obs.qcindex  ) = obs.qc;
obsmat(:,obs.keyindex ) = obs.keys;
obsmat(:,obs.timeindex) = obs.time;
obsmat(:,obs.indindex ) = [1:length(obs.time)];

linked_observations(obsmat,obs)
