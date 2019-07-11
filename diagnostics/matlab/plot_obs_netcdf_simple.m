function obsstruct = plot_obs_netcdf_simple(fname, varargin)

%% plot_obs_netcdf_simple(fname) will plot the locations and values of the observations in a DART netcdf file.
%     There are several [PARAM,VALUE] optional arguments for more control over what gets plotted.
%
% plot_obs_netcdf_simple(fname,'region',[180 200 0 90 -Inf Inf]); % will restrict the region of what is plotted.
%
% a list of PARAM,VALUE pairs and their defaults follows:
%
% ObsTypeString = 'ALL';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'NCEP BUFR observation';
% QCString      = 'DART quality control';
% maxgoodQC     = 0;   % in this context, anything with a DART QC > 0 is plotted in a separate axes
% verbose       = 0;   % anything > 0 == 'true'
% twoup         = 0;   % anything > 0 means plot both axes on same figure
%
% todo ... if there is only one copy or one qc string - just use it regardless of what has been input

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Provide defaults, Decode, Parse, Check the input
%---------------------------------------------------------------------

default_ObsTypeString = 'ALL';
default_region        = [0 360 -90 90 -Inf Inf];
default_CopyString    = 'NCEP BUFR observation';
default_QCString      = 'DART quality control';
default_maxgoodQC     = 0;
default_verbose       = 0;
default_twoup         = 0;
p = inputParser;

addRequired(p,'fname',@ischar);
% addParameter or addParamValue depend on version of matlab - support both
if (exist('inputParser/addParameter','file') == 2)
    addParameter( p,'ObsTypeString',default_ObsTypeString,@ischar);
    addParameter( p,'region',       default_region       ,@isnumeric);
    addParameter( p,'CopyString',   default_CopyString   ,@ischar);
    addParameter( p,'QCString',     default_QCString     ,@ischar);
    addParameter( p,'maxgoodQC',    default_maxgoodQC    ,@isnumeric);
    addParameter( p,'verbose',      default_verbose      ,@isnumeric);
    addParameter( p,'twoup',        default_twoup        ,@isnumeric);
else
    addParamValue(p,'ObsTypeString',default_ObsTypeString,@ischar);
    addParamValue(p,'region',       default_region       ,@isnumeric);
    addParamValue(p,'CopyString',   default_CopyString   ,@ischar);
    addParamValue(p,'QCString',     default_QCString     ,@ischar);
    addParamValue(p,'maxgoodQC',    default_maxgoodQC    ,@isnumeric);
    addParamValue(p,'verbose',      default_verbose      ,@isnumeric);
    addParamValue(p,'twoup',        default_twoup        ,@isnumeric);
end

p.parse(fname, varargin{:});

% if you want to echo the input
% disp(['fname         : ', p.Results.fname])
% disp(['p.Results.ObsTypeString : ', p.Results.p.Results.ObsTypeString])

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

if (exist(fname,'file') ~= 2)
    error('%s does not exist.',fname)
end

if ( p.Results.twoup > 0 )
    clf; orient tall
    positions = [0.1, 0.55, 0.8, 0.35 ; ...
        0.1, 0.10, 0.8, 0.35 ; ...
        0.1, 0.02, 0.8, 0.08];
else
    clf; orient landscape
    positions = [0.1, 0.20, 0.8, 0.65 ; ...
        0.1, 0.20, 0.8, 0.65 ; ...
        0.1, 0.05, 0.8, 0.10];
end

%% Read the observation sequence

obsstruct = read_obs_netcdf(p.Results.fname, p.Results.ObsTypeString, p.Results.region, ...
    p.Results.CopyString, p.Results.QCString, p.Results.verbose);

% subset based on qc value

if ( (~ isempty(obsstruct.qc)) && (~ isempty(p.Results.maxgoodQC)) )

    inds = find(obsstruct.qc > p.Results.maxgoodQC);
    obsstruct.numflagged = length(inds);
    fprintf('%d obs have a %s value greater than %f (i.e. "bad")\n', ...
        obsstruct.numflagged, p.Results.QCString, p.Results.maxgoodQC)

    if (~isempty(inds))
        flaggedobs.lons = obsstruct.lons(inds);
        flaggedobs.lats = obsstruct.lats(inds);
        flaggedobs.Ztyp = obsstruct.Ztyp(inds);
        flaggedobs.z    = obsstruct.z(   inds);
        flaggedobs.obs  = obsstruct.obs( inds);
        flaggedobs.qc   = obsstruct.qc(  inds);
    end

    inds = find(obsstruct.qc <= p.Results.maxgoodQC);
    obsstruct.numgood = length(inds);
    fprintf('%d obs have a %s value less than or equal to %f (i.e. "good")\n', ...
        obsstruct.numgood, p.Results.QCString, p.Results.maxgoodQC)

    bob = obsstruct.lons(inds); obsstruct.lons = bob;
    bob = obsstruct.lats(inds); obsstruct.lats = bob;
    bob = obsstruct.Ztyp(inds); obsstruct.Ztyp = bob;
    bob = obsstruct.z(   inds); obsstruct.z    = bob;
    bob = obsstruct.obs( inds); obsstruct.obs  = bob;
    bob = obsstruct.qc(  inds); obsstruct.qc   = bob;

end

%% Separate out all the NaN values
%  If there are NaN's that make it through to the plotting phase,
%  the entire plot gets wiped out when the continents get plotted.
%  The locations with NaN values will be plotted with separate symbols.

numobs = length(obsstruct.obs);
ngood  = sum(isfinite(obsstruct.obs));

if numobs == ngood

    nanobs.numflagged = 0;
    nanobs.string     = 'No "good" observations with "NaN" values.';

else

    % Creating a structure with the unplottable values.

    inds = find(isfinite(obsstruct.obs) == 0);

    nanobs.numflagged = length(inds);
    nanobs.string = sprintf('%d obs have unplottable values.\n', nanobs.numflagged);

    if (~isempty(inds))
        nanobs.lons = obsstruct.lons(inds);
        nanobs.lats = obsstruct.lats(inds);
        nanobs.Ztyp = obsstruct.Ztyp(inds);
        nanobs.z    = obsstruct.z(   inds);
        nanobs.obs  = 1;
        nanobs.qc   = obsstruct.qc(  inds);
    end

    % Removing the unplottable values from the plottable ones.
    fprintf('%d obs have unplottable values.\n', nanobs.numflagged)

    inds = isfinite(obsstruct.obs);

    bob = obsstruct.lons(inds); obsstruct.lons = bob;
    bob = obsstruct.lats(inds); obsstruct.lats = bob;
    bob = obsstruct.Ztyp(inds); obsstruct.Ztyp = bob;
    bob = obsstruct.z(   inds); obsstruct.z    = bob;
    bob = obsstruct.obs( inds); obsstruct.obs  = bob;
    bob = obsstruct.qc(  inds); obsstruct.qc   = bob;

end

%% Create graphic with area-weighted symbols for the good observations.
%  It has happened that there have been zero good observations in a file.

xmin = min(p.Results.region(1:2));
xmax = max(p.Results.region(1:2));
ymin = min(p.Results.region(3:4));
ymax = max(p.Results.region(3:4));
zmin = min(obsstruct.z);
zmax = max(obsstruct.z);

pstruct.colorbarstring = obsstruct.ObsTypeString;
pstruct.p.Results.region = p.Results.region;
pstruct.str1   = sprintf('%s',obsstruct.ObsTypeString);
pstruct.str3   = sprintf('%s - %s',obsstruct.timestring(1,:),obsstruct.timestring(2,:));

subplot('position',positions(1,:))

if ( length(obsstruct.obs) < 1 )
    fprintf('There are no ''good'' observations to plot.\n')
    % may still want to plot the obs with bad qc values.
    str1 = sprintf('There are no ''good'' observations to plot.\n');
    text(0.5,0.67,str1,'HorizontalAlignment','center')
    text(0.5,0.33,nanobs.string,'HorizontalAlignment','center')
    title( {pstruct.str1, obsstruct.CopyString, pstruct.str3}, 'Interpreter','none','FontSize',14);
else

    pstruct.Ztype  = obsstruct.Ztyp(1);

    % choose a symbol size based on the number of obs to plot.
    if (length(obsstruct.obs) < 1000)
        pstruct.scalearray = scaleme(obsstruct.obs, 30);
    else
        pstruct.scalearray = 30.0 * ones(size(obsstruct.obs));
    end
    pstruct.clim   = [min(obsstruct.obs)-0.5 max(obsstruct.obs)+0.5];
    pstruct.str2   = sprintf('%s (%d locations)',obsstruct.CopyString,length(obsstruct.obs));

    % If all the observations live on the same level ... make a 2D plot.

    if ( zmin == zmax )

        pstruct.axis = [xmin xmax ymin ymax zmin-0.5 zmax+0.5];

        plot_2D(obsstruct, pstruct);

    else

        pstruct.axis = [xmin xmax ymin ymax zmin zmax];
        pstruct.str1 = sprintf('%s level (%.2f - %.2f)',obsstruct.ObsTypeString,zmin,zmax);

        plot_3D(obsstruct, pstruct);

    end

    % If there are unplottable values ... report that.

    if nanobs.numflagged > 0
        subplot('position',positions(3,:))
        axis off
        text(0.5,0.5,nanobs.string,'HorizontalAlignment','center')
    end

end

%% Create graphic of spatial distribution of 'flagged' observations & their QC value.
%
% 0     observation assimilated
% 1     observation evaluated only
%   --- everything above this means the prior and posterior are OK
% 2     assimilated, but the posterior forward operator failed
% 3     Evaluated only, but the posterior forward operator failed
%   --- everything above this means only the prior is OK
% 4     prior forward operator failed
% 5     not used
% 6     prior QC rejected
% 7     outlier rejected

dartqc_strings = { ...
    '''observation evaluated only''', ...
    '''assimilated, but the posterior forward operator failed''', ...
    '''evaluated only, but the posterior forward operator failed''',...
    '''prior forward operator failed''',...
    '''not used''',...
    '''prior QC rejected''',...
    '''outlier rejected''',...
    '''vertical conversion failed''', ...
    '''reserved for future use'''};

if (obsstruct.numflagged > 0 ) % if there are flagged observation to plot ... carry on.

    pstruct.Ztype  = flaggedobs.Ztyp(1);

    if (p.Results.twoup <= 0)
        figure; clf
    end

    subplot('position',positions(2,:))

    zmin = min(flaggedobs.z);
    zmax = max(flaggedobs.z);

    prej = 100.0 * length(flaggedobs.obs) / ...
        (length(flaggedobs.obs) + length(obsstruct.obs));

    % choose a symbol size based on the number of obs to plot.
    if (length(flaggedobs.obs) < 1000)
        pstruct.scalearray = scaleme(flaggedobs.obs, 30);
    else
        pstruct.scalearray = 30.0 * ones(size(flaggedobs.obs));
    end

    pstruct.colorbarstring = p.Results.QCString;
    pstruct.clim = [min(flaggedobs.qc) max(flaggedobs.qc)];
    pstruct.str1 = sprintf('%s level (%.2f - %.2f)',obsstruct.ObsTypeString,zmin,zmax);
    pstruct.str2 = sprintf('%s (%d ''good'', %d ''flagged'' -- %.2f %%)', obsstruct.CopyString, ...
        length(obsstruct.obs), length(flaggedobs.obs), prej);

    flaggedobs.obs = flaggedobs.qc;  % plot QC values, not obs values
    if ( zmin ~= zmax )

        pstruct.axis = [xmin xmax ymin ymax zmin zmax];

        plot_3D(flaggedobs, pstruct);

    else

        pstruct.axis = [xmin xmax ymin ymax zmin-0.5 zmin+0.5];

        plot_2D(flaggedobs, pstruct);

    end

    subplot('position',positions(3,:))
    axis off

    %% If the QC is a DART QC, we know how to interpret them.

    switch lower(strtrim(p.Results.QCString))
        case 'dart quality control',

            qcvals  = unique(flaggedobs.qc);
            qccount = zeros(size(qcvals));
            s = cell(length(qcvals));
            for i = 1:length(qcvals)
                qccount(i) = sum(flaggedobs.qc == qcvals(i));
                s{i} = sprintf('%d obs with qc == %d %s',qccount(i),qcvals(i), ...
                    dartqc_strings{qcvals(i)});
            end

            dy =  0.8*1.0/length(s);
            for i = 1:length(s)
                text(0.0, (i-1)*dy ,s{i})
            end

        otherwise,
            str = sprintf('no way to interpret values of %s',strtrim(p.Results.QCString));
            text(0.0, 0.0, str)
    end
end



function s = scaleme(x,minsize)
% scaleme returns a uniformly scaled array the same size as the input
% array where the maximum is 10 times the minimum
maxsize = 10*minsize;
minx    = min(x);
maxx    = max(x);
slope   = (maxsize-minsize)/(maxx-minx);
b       = minsize - slope*minx;

s = x*slope + b;



function h1 = plot_3D(obsstruct, pstruct)

scatter3(obsstruct.lons, obsstruct.lats, obsstruct.z, ...
    pstruct.scalearray, obsstruct.obs, 'd', 'filled');

h1 = gca;

axis(pstruct.axis)

title( {pstruct.str1, pstruct.str3, pstruct.str2}, 'Interpreter','none','FontSize',14);
xlabel('longitude')
ylabel('latitude')

if     (obsstruct.Ztyp(1) == -2) % VERTISUNDEF     = -2
   zlabel('unspecified')
elseif (obsstruct.Ztyp(1) == -1) % VERTISSURFACE   = -1
   zlabel('surface')
elseif (obsstruct.Ztyp(1) ==  1) % VERTISLEVEL     =  1
   zlabel('level')
elseif (obsstruct.Ztyp(1) ==  2) % VERTISPRESSURE  =  2
   set(gca,'ZDir','reverse')
   zlabel('pressure')
elseif (obsstruct.Ztyp(1) ==  3) % VERTISHEIGHT    =  3
   zlabel('height')
end

FlatEarth(pstruct);
hb = colorbar;
set(get(hb,'YLabel'),'String',pstruct.colorbarstring,'Interpreter','none')




function h1 = plot_2D(obsstruct, pstruct)

if (isfield(obsstruct,'z'))
    scatter3(obsstruct.lons, obsstruct.lats, obsstruct.z, ...
        pstruct.scalearray, obsstruct.obs, 'd', 'filled');
else
    scatter3(obsstruct.lons, obsstruct.lats, zeros(size(obsstruct.lons)), ...
        pstruct.scalearray, obsstruct.obs, 'd', 'filled');
end

h1 = gca;

axis(pstruct.axis);

title( {pstruct.str1, pstruct.str3, pstruct.str2}, 'Interpreter','none','FontSize',14);
xlabel('longitude')
ylabel('latitude')

FlatEarth(pstruct);
hb = colorbar;
set(get(hb,'YLabel'),'String',pstruct.colorbarstring,'Interpreter','none')
view(0,90)


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
