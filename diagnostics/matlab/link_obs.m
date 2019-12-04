function link_obs(fname, ObsTypeString, ObsCopyString, CopyString, QCString, region)
%% link_obs generates the 'brushable' observation plots.
%
% link_obs explores the observations in a netCDF file created by 'obs_seq_to_netcdf'.
% All the data in the graphics are linked and are 'brushable'. Click on
% the paintbrush icon and use the mouse to select observations in any
% graphic. That same observation will be highlighted in ALL the graphics.
%
% Three figures will be generated.
%
% Figure 1 will have a 3D geographic scatterplot. Click on the rotate
%          icon and drag the graphic around for the best view angle.
%
% Figure 2 has multiple plots.
%   *   The bottom plot is the observation 'key'
%       (some of the linked list information in the original file).
%   *   The middle plot provides information about the observation
%       density as a function of time.
%   *   The final (top) plot has the QC value as a function of time.
%
% Figure 3 has two plots.
%   *   The bottom plot is a 2D scatterplot of (typically) the
%       'prior mean observation' (or whatever you specify in Copystring) vs. the
%       'observation' (or whatever you specify in ObsCopyString).
%       Both of these can be changed however - the allowable settings are defined
%       in the CopyMetaData variable in the netCDF file.
%
% link_obs(fname, ObsTypeString, ObsCopyString, CopyString, QCString, region)
%
% fname         - name of netcdf file that results from obs_seq_to_netcdf
% ObsTypeString - the TYPE of observation (ncdump -v ObsTypesMetaData *.nc)
% ObsCopyString - the COPY specifying the raw observation ( -v CopyMetaData )
% CopyString    - the COPY specifying the copy to compare to the raw obs
% QCString      - character string  - one of (ncdump -v QCMetaData *.nc)
% region        - geographic extent of interest
%
% EXAMPLE 1:
% fname         = '/ptmp/thoar/POP/CAM/POP8/obs_epoch_001.nc';
% ObsTypeString = 'APB_TEMPERATURE';
% ObsCopyString = 'WOD observation';
% CopyString    = 'prior ensemble mean';
% QCString      = 'DART quality control';
% region        = [0 360 -90 90 -Inf Inf];
% global obsmat;
% link_obs(fname, ObsTypeString, ObsCopyString, CopyString, QCString, region)
%
% EXAMPLE 2:
% fname         = 'obs_epoch_001.nc';
% ObsTypeString = 'RADIOSONDE_TEMPERATURE';
% ObsCopyString = 'NCEP BUFR observation';
% CopyString    = 'prior ensemble mean';
% QCString      = 'DART quality control';
% region        = [220 300 20 60 -Inf Inf];
% global obsmat;
% link_obs(fname, ObsTypeString, ObsCopyString, CopyString, QCString, region)
%
% EXAMPLE 3:
% global obsmat;
% region        = [0 360 -90 90 -Inf Inf];
% link_obs('obs_epoch_002.nc','RADIOSONDE_TEMPERATURE', 'observation', ...
%          'prior ensemble member 3', 'DART quality control', region)
%
% IMPORTANT: click on the little paintbrush icon in order to activate the
% 'brushable' feature on the plots. Once that is highlighted, any observations
% selected in one view become highlighted in ALL the views.
% REALLY IMPORTANT: pre-2018b  the paintbrush icon is on the figure toolbar.
% REALLY IMPORTANT: 2018b-?  the paintbrush icon is on the AXES toolbar, 
%         you have to hover over a plot and the axes toolbar 'appears'.
%
% ALSO IMPORTANT: If you are using the Matlab GUI, doubleclick on 'obsmat' in
% the the Workspace window to generate a spreadsheet-like view of the
% observations which is also linked to the data brushing.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (exist(fname,'file') ~= 2)
    error('%s does not exist.',fname)
end

verbose = 1;

obs  = read_obs_netcdf(fname, ObsTypeString, region, ObsCopyString, QCString, verbose);
copy = read_obs_netcdf(fname, ObsTypeString, region,    CopyString, QCString, 0);

if ( isempty(obs.lons) )
    error('There are no %s observations in the region specified in %s', ObsTypeString, fname)
end

l1 = length(copy.obs);
l2 = length(obs.obs);

if (l1 ~= l2 )
    error('there are %d %s and only %d %s',l1, CopyString, l2, ObsCopyString)
end

obs.ObsCopyString = obs.CopyString;
obs.CopyString    = copy.CopyString;
obs.region(5)     = min(obs.z); % use observation Z to specify vertical region
obs.region(6)     = max(obs.z);

% Some observations are all at same z value. Cannot use same value for
% upper and lower axis bounds.
if (obs.region(5) == obs.region(6))
    obs.region(5) = obs.region(5) - 0.5;
    obs.region(6) = obs.region(6) + 0.5;
end

%% Now pack the data in the same fashion as the cell array of column labels.

obs.lonindex  = 1;
obs.latindex  = 2;
obs.zindex    = 3;
obs.obsindex  = 4;
obs.copyindex = 5;
obs.qcindex   = 6;
obs.keyindex  = 7;
obs.timeindex = 8;
obs.indindex  = 9;

global obsmat
obsmat = zeros(length(obs.lons),9);
obsmat(:,obs.lonindex ) = obs.lons; obs.colnames{obs.lonindex}  = 'longitude';
obsmat(:,obs.latindex ) = obs.lats; obs.colnames{obs.latindex}  = 'latitude';
obsmat(:,obs.zindex   ) = obs.z   ; obs.colnames{obs.zindex}    = obs.Zunits;
obsmat(:,obs.obsindex ) = obs.obs ; obs.colnames{obs.obsindex}  = ObsCopyString;
obsmat(:,obs.copyindex) = copy.obs; obs.colnames{obs.copyindex} = CopyString;
obsmat(:,obs.qcindex  ) = obs.qc  ; obs.colnames{obs.qcindex}   = QCString;
obsmat(:,obs.keyindex ) = obs.keys; obs.colnames{obs.keyindex}  = 'obskey';
obsmat(:,obs.timeindex) = obs.time; obs.colnames{obs.timeindex} = 'time';
obsmat(:,obs.indindex ) = 1:length(obs.time); obs.colnames{obs.indindex} = 'index';

%% Replace all ill-posed copies with
% This should really check to see if the copy is a 'mean' or 'spread' ...

iscalculated = ~isempty(strfind(lower(obs.colnames{obs.copyindex}),'mean')) | ...
    ~isempty(strfind(lower(obs.colnames{obs.copyindex}),'spread'));

if iscalculated
    disp('replacing copies with [1 < QC flag < 5] with NaN')
    allindices = obsmat(:,obs.qcindex);
    badinds = ((allindices > 1) & (allindices < 5));
    obsmat(badinds,obs.copyindex) = NaN;
end

%% create the linked plots
linked_observations(obs)

fprintf('\nUse the paintbrush icon to select data.\n')
if verLessThan('matlab','R2018b') 
   fprintf('The paintbrush is on any figure toolbar.\n\n')
else
   fprintf('The paintbrush is available by hovering over a plot.\n')
   fprintf('It will turn blue when enabled.\n\n')
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
