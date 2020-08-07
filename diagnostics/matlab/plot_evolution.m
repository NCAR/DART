function plotdat = plot_evolution(fname, copy, varargin)
%% plot_evolution plots the temporal evolution of the observation-space quantities for all possible levels, all possible variables.
% Part of the observation-space diagnostics routines.
%
% 'obs_diag' produces a netcdf file containing the diagnostics.
% 'obs_diag' condenses the obs_seq.final information into summaries for a few
% specified regions - on a level-by-level basis.
%
% The number of observations possible reflects only those observations
% that have incoming QC values of interest. Any observation with a DART
% QC of 5 or 6 is not considered 'possible' for the purpose of this graphic.
%
% NOTE: if the observation was designated as a TRUSTED observation in the
%       obs_diag program, the observations that were rejected by the outlier
%       threshhold STILL PARTICIPATE in the calculation of the rmse, spread, etc.
%       The _values_ plotted by plot_profile reflect that. The number of observations
%       "used" becomes unclear. The number of observations used (designated by the
%       asterisk) is ALWAYS the number of observations successfully assimilated.
%       For TRUSTED observations, this is different than the number used to calculate
%       bias, rmse, spread, etc.
%
% USAGE: plotdat = plot_evolution(fname, copy [,varargin]);
%
% fname    :  netcdf file produced by 'obs_diag'
%
% copy     : string defining the metric of interest. 'rmse', 'spread', etc.
%            Possible values are available in the netcdf 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% varargin: optional parameter-value pairs. Supported parameters are described below.
%
% obsname  : The strings of each observation type to plot.
%            Each observation type will be plotted in a separate graphic.
%            Default is to plot all available observation types.
%
% level    : 'level' index. Default is to plot all levels.
%
% range    : 'range' of the value being plotted. Default is to
%                automatically determine range based on the data values.
%
% verbose  : true/false to control amount of run-time output
%
% MarkerSize : integer controlling the size of the symbols
%
% DateForm : Free-form character string controlling representation of the time axis.
%            See 'help datetick' for discussion and valid values.
%            Example ones are 'mm/dd' and 'dd HH:MM'.
%
% pause  : true/false to conrol pausing after each figure is created.
%          true will require hitting any key to continue to next plot
%
% OUTPUT: 'plotdat' is a structure containing what was last plotted.
%         A postscript file containing a page for each level - each region.
%         The other file is a simple text file containing summary information
%         about how many observations were assimilated, how many were available, etc.
%         Both of these filenames contain the observation type,
%         copy and region as part of the name.
%
%
% EXAMPLE 1 - plot the evolution of the bias for all observation types, all levels
%
% fname   = 'obs_diag_output.nc';
% copy    = 'bias';
% plotdat = plot_evolution(fname, copy);
%
%
% EXAMPLE 2 - plot the evolution of the rmse for just the radiosonde temperature obs
%             This requires that the 'RADIOSONDE_TEMPERATURE' is one of the known
%             observation types in the netCDF file.
%
% fname   = 'obs_diag_output.nc';
% copy    = 'rmse';
% plotdat = plot_evolution(fname, copy, 'obsname', 'RADIOSONDE_TEMPERATURE');
%
%
% EXAMPLE 3 - plot the evolution of the rmse for just the radiosonde temperature obs
%             for the 4th level and force the vertical axis of the 'rmse' to be 0,10
%
% plotdat    = plot_evolution(fname, 'rmse', 'obsname', 'RADIOSONDE_TEMPERATURE', ...
%                             'level', 4, 'range', [0 10], 'pause', false);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

default_obsname    = 'none';
default_verbosity  = true;
default_markersize = 12;
default_pause      = false;
default_range      = [NaN NaN];
default_level      = -1;
default_dateform   = 'default';
p = inputParser;

addRequired(p,'fname',@ischar);
addRequired(p,'copy',@ischar);
if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'obsname',    default_obsname,    @ischar);
    addParameter(p,'verbose',    default_verbosity,  @islogical);
    addParameter(p,'MarkerSize', default_markersize, @isnumeric);
    addParameter(p,'pause',      default_pause,      @islogical);
    addParameter(p,'range',      default_range,      @isnumeric);
    addParameter(p,'level',      default_level,      @isnumeric);
    addParameter(p,'DateForm',   default_dateform,   @ischar);
else
    addParamValue(p,'obsname',   default_obsname,    @ischar);    %#ok<NVREPL>
    addParamValue(p,'verbose',   default_verbosity,  @islogical); %#ok<NVREPL>
    addParamValue(p,'MarkerSize',default_markersize, @isnumeric); %#ok<NVREPL>
    addParamValue(p,'pause',     default_pause,      @islogical); %#ok<NVREPL>
    addParamValue(p,'range',     default_range,      @isnumeric); %#ok<NVREPL>
    addParamValue(p,'level',     default_level,      @isnumeric); %#ok<NVREPL>
    addParamValue(p,'DateForm',  default_dateform,   @ischar);    %#ok<NVREPL>
end
p.parse(fname, copy, varargin{:});

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

if (numel(p.Results.range) ~= 2)
    error('range must be an array of length two ... [bottom top]')
end

if strcmp(p.Results.obsname,'none')
    nvars = 0;
else
    obsname = p.Results.obsname;
    nvars = 1;
end

if (exist(fname,'file') ~= 2)
    error('file/fname <%s> does not exist',fname)
end

%%--------------------------------------------------------------------
% Harvest plotting info/metadata from netcdf file.
%---------------------------------------------------------------------

plotdat = read_obsdiag_staticdata(fname,copy);

if (nvars == 0)
    [plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(fname);
    [plotdat.varnames,    plotdat.vardims]    = FindTemporalVars(plotdat);
    plotdat.nvars       = length(plotdat.varnames);
else
    plotdat.varnames{1} = obsname;
    plotdat.nvars       = nvars;
end

global figuredata verbose

figuredata            = set_obsdiag_figure('landscape');
figuredata.MarkerSize = p.Results.MarkerSize;
figuredata.DateForm   = p.Results.DateForm;
verbose               = p.Results.verbose;

%%---------------------------------------------------------------------
% Loop around (time-copy-level-region) observation types
%----------------------------------------------------------------------

psfname = cell(plotdat.nvars);

for ivar = 1:plotdat.nvars
    
    % create the variable names of interest.
    
    plotdat.myvarname = plotdat.varnames{ivar};
    plotdat.guessvar  = sprintf('%s_guess',plotdat.varnames{ivar});
    plotdat.analyvar  = sprintf('%s_analy',plotdat.varnames{ivar});
    
    plotdat.trusted   = nc_read_att(fname, plotdat.guessvar, 'TRUSTED');
    if (isempty(plotdat.trusted)), plotdat.trusted = 'NO'; end
    
    % remove any existing postscript file - will simply append each
    % level as another 'page' in the .ps file.
    
    for iregion = 1:plotdat.nregions
        psfname{iregion} = sprintf('%s_%s_evolution_region%d.ps', ...
            plotdat.varnames{ivar}, plotdat.copystring, iregion);
        if (exist(psfname{iregion},'file') == 2)
            fprintf('Removing %s from the current directory.\n',psfname{iregion})
            system(sprintf('rm %s',psfname{iregion}));
        end
    end
    
    % remove any existing log file -
    
    lgfname = sprintf('%s_%s_obscount.txt',plotdat.varnames{ivar},plotdat.copystring);
    if (exist(lgfname,'file') == 2)
        fprintf('Removing %s from the current directory.\n',lgfname)
        system(sprintf('rm %s',lgfname));
    end
    logfid = fopen(lgfname,'wt');
    fprintf(logfid,'%s\n',lgfname);
    
    % check to see if there is anything to plot
    % The number possible is decreased by the number of observations
    % rejected by namelist control.
    
    qcvalues = get_qc_values(fname, plotdat.guessvar, 'fatal', false, ...
        'verbose', false);
    
    if ( sum(qcvalues.nposs(:)) < 1 )
        fprintf('no obs for %s...  skipping\n', plotdat.varnames{ivar})
        continue
    end
    
    % get appropriate vertical coordinate variable
    
    [dimnames, ~] = nc_var_dims(fname, plotdat.guessvar);
    
    if ( plotdat.dimensionality == 1 ) % observations on a unit circle, no level
        plotdat.varlevels   = 1;
        plotdat.level       = 1;
        plotdat.level_units = [];
    elseif ( strfind(dimnames{2},'surface') > 0 )
        plotdat.varlevels   = 1;
        plotdat.level       = 1;
        plotdat.level_units = 'surface';
    elseif ( strfind(dimnames{2},'undef') > 0 )
        plotdat.varlevels   = 1;
        plotdat.level       = 1;
        plotdat.level_units = 'undefined';
    else
        plotdat.varlevels   = ncread(fname, dimnames{2});
        plotdat.level_units = nc_read_att(fname, dimnames{2}, 'units');
        nlevels             = length(plotdat.varlevels);
        if (p.Results.level < 0 )
            % use all the levels
            plotdat.level   = 1:nlevels;
        elseif (p.Results.level > 0 && p.Results.level < nlevels)
            plotdat.level   = round(p.Results.level);
        else
            str1 = sprintf('valid level values are 1 <= %d',nlevels);
            error('%s\n%f is not a valid level for %s\n',str1,p.Results.level,plotdat.guessvar)
        end
    end
    
    % read the whole variable, subset it later
    
    guess = ncread(fname, plotdat.guessvar);
    analy = local_ncread(fname, plotdat.analyvar);
    if ( isempty(analy) )
        % force analysis to be the same shape as the guess and full of NaNs
        analy = guess;
        analy(:) = NaN;
        has_posterior = false;
        plotdat.post_string = '';
    else
        has_posterior = true;
        plotdat.post_string = '; \diamondsuit=posteriorOK';
    end
    
    for ilevel = 1:length(plotdat.level)
        priorQCs = get_qc_values(fname, plotdat.guessvar, ...
            'levelindex', plotdat.level(ilevel), ...
            'fatal', false, ...
            'verbose', verbose);
        plotdat.mylevel   = plotdat.level(ilevel);
        plotdat.ges_Neval = priorQCs.num_evaluated;
        plotdat.ges_Nposs = priorQCs.nposs;
        plotdat.ges_Nused = priorQCs.nused;
        plotdat.ges_copy  = guess(:,plotdat.mylevel,plotdat.copyindex,:);
        plotdat.anl_copy  = analy(:,plotdat.mylevel,plotdat.copyindex,:);
        
        if (has_posterior)
            posteQCs = get_qc_values(fname, plotdat.analyvar, ...
                'levelindex', plotdat.mylevel, ...
                'fatal', false, ...
                'verbose', verbose);
            plotdat.anl_Nused = posteQCs.nused;
            plotdat.anl_copy  = analy(:,plotdat.mylevel,plotdat.copyindex,:);
        else
            plotdat.anl_Nused = zeros(size(plotdat.ges_Nused));
            plotdat.anl_copy  = plotdat.ges_copy;  % needed for determining limits
        end
        
        % call report_qc_values.m
        
        if isnan(p.Results.range(1))
            plotdat.Yrange = FindRange(plotdat);
        else
            plotdat.Yrange = p.Results.range;
        end
        
        % plot each region, each level to a separate figure
        
        for iregion = 1:plotdat.nregions
            figure(iregion); clf(iregion); orient(figuredata.orientation);
            
            plotdat.region   = iregion;
            plotdat.myregion = deblank(plotdat.region_names(iregion,:));
            if ( isempty(plotdat.level_units) )
                plotdat.title = plotdat.myvarname;
            else
                plotdat.title = sprintf('%s @ %f %s',    ...
                    plotdat.myvarname,     ... 
                    plotdat.varlevels(plotdat.mylevel), ...
                    plotdat.level_units);
            end
            
            myplot(plotdat);
            
            % create/append to the postscript file
            if verLessThan('matlab','R2016a')
                print(gcf, '-dpsc', '-append', psfname{iregion});
            else
                print(gcf, '-dpsc', '-append', '-bestfit', psfname{iregion});
            end
            
            % block to go slow and look at each one ...
            if (p.Results.pause)
                disp('Pausing, hit any key to continue ...')
                pause
            end
            
        end
    end
end

%=====================================================================
% 'Helper' functions
%=====================================================================


function myplot(plotdat)

%% The prior and posterior are plotted as separate items.
% By this point, the middle two dimensions are singletons.

global figuredata verbose

ax1 = subplot('position',figuredata.position);
set(ax1,'YAxisLocation','left','FontSize',figuredata.fontsize)
orient(figuredata.orientation)

[hprior, prior_legstr] = plot_quantity('prior', plotdat);

ges_Nposs = squeeze(plotdat.ges_Nposs(plotdat.region,:,:,:));
ges_Nused = squeeze(plotdat.ges_Nused(plotdat.region,:,:,:));
anl_Nused = squeeze(plotdat.anl_Nused(plotdat.region,:,:,:));
anl_Ngood = sum(anl_Nused);

if anl_Ngood
    [hposte, poste_legstr] = plot_quantity('posterior', plotdat);
    h = legend([hprior, hposte], prior_legstr, poste_legstr);
else
    h = legend(hprior,prior_legstr);
    poste_legstr = [];
end

set(h,'Interpreter','none','Box','off','FontSize',figuredata.fontsize)

if verLessThan('matlab','R2017a')
    % Convince Matlab to not autoupdate the legend with each new line.
    % Before 2017a, this was the default behavior, so do nothing.
    % We do not want to add the bias line to the legend, for example.
else
    h.AutoUpdate = 'off';
end

if verbose
    fprintf('region %d %s level %f nobs_poss %d prior %d poste %d\n', ...
        plotdat.region, plotdat.myvarname, plotdat.mylevel, ...
        sum(ges_Nposs), sum(ges_Nused), anl_Ngood)
    fprintf('region %d %s level %f %s %s\n\n', ...
        plotdat.region, plotdat.myvarname, plotdat.mylevel, prior_legstr, poste_legstr)
end

% Attempt to make plotting robust in the face of 'empty' bins.
% The bincenters variable has all the temporal bins specified,
% so we use that to determine the X axis limits.

axlims = [min(plotdat.bincenters) max(plotdat.bincenters) plotdat.Yrange];
axis(axlims)

switch lower(plotdat.copystring)
    case 'bias'
        % plot a zero-bias line
        zeroline = line(axlims(1:2),[0 0], 'Color',[200 200 200]/255,'Parent',ax1);
        set(zeroline,'LineWidth',2.5,'LineStyle','-')
        plotdat.ylabel = sprintf('%s (%s)',plotdat.copystring,plotdat.biasconv);
    otherwise
        plotdat.ylabel = sprintf('%s',plotdat.copystring);
end

% effort to use user-supplied value for time labelling or
% make a stab at a useful default.

xlabelstring = set_time_axis('x', plotdat.bincenters, figuredata.DateForm);

set(get(ax1,'Xlabel'),'String',xlabelstring, ...
    'Interpreter','none','FontSize',figuredata.fontsize)

title({plotdat.myregion, plotdat.title}, ...
    'Interpreter', 'none', 'Fontsize', figuredata.fontsize, 'FontWeight', 'bold')
BottomAnnotation(plotdat)

% create a separate scale for the number of observations
ax2 = axes( ...
    'Position',get(ax1,'Position'), ...
    'FontSize',get(ax1,'FontSize'), ...
    'XColor'  ,get(ax1,'XColor'), ...
    'XLim'    ,get(ax1,'XLim'), ...
    'XTick'   ,get(ax1,'XTick'), ...
    'YDir'    ,get(ax1,'YDir'), ...
    'Color'   ,'none', ...
    'YColor'  ,figuredata.obs_color, ...
    'XAxisLocation','top', ...
    'YAxisLocation','right');

ax2h1 = line(plotdat.bincenters, ges_Nposs, 'Parent', ax2);
ax2h2 = line(plotdat.bincenters, ges_Nused, 'Parent', ax2);

set(ax2h1, 'LineStyle', 'none', ...
    'Color',     figuredata.obs_color, ...
    'Marker',    figuredata.obs_marker, ...
    'MarkerSize',figuredata.MarkerSize);

set(ax2h2, 'LineStyle', 'none', ...
    'Color',     figuredata.obs_color, ...
    'Marker',    figuredata.ges_marker, ...
    'MarkerSize',figuredata.MarkerSize);

if anl_Ngood > 0
    ax2h3 = line(plotdat.bincenters, anl_Nused, 'Parent',ax2);
    set(ax2h3, 'LineStyle', 'none', ...
        'Color',     figuredata.obs_color, ...
        'Marker',    figuredata.anl_marker, ...
        'MarkerSize',figuredata.MarkerSize);
end

% turn off topside X tick labels (clashes with title)
% use the same Y ticks, but find the right label values
set(ax2, 'XTicklabel', []);
matchingYticks(ax1,ax2);

set(get(ax1,'Ylabel'), 'String', plotdat.ylabel, ...
    'Interpreter','none','FontSize',figuredata.fontsize)

% determine if the observation type was flagged as 'evaluate' or 'assimilate'
% since we don't have the ability to specify this level-by-level or by
% regions, we can use an 'all-or-nothing' approach.

if sum(plotdat.ges_Neval(:)) > 0
    string1 = ['# of obs: o=possible; \ast=evaluated' plotdat.post_string];
else
    string1 = ['# of obs: o=possible; \ast=assimilated' plotdat.post_string];
end
set(get(ax2,'Ylabel'), 'String', string1, 'FontSize', figuredata.fontsize)


%=====================================================================


function BottomAnnotation(main)
%% annotates the full path of the file being plotted
subplot('position',[0.48 0.01 0.04 0.04])
axis off
fullname = which(main.fname);   % Could be in MatlabPath
if( isempty(fullname) )
    if ( main.fname(1) == '/' )  % must be a absolute pathname
        string1 = sprintf('data file: %s',main.fname);
    else                   % must be a relative pathname
        mydir = pwd;
        string1 = sprintf('data file: %s/%s',mydir,main.fname);
    end
else
    string1 = sprintf('data file: %s',fullname);
end

h = text(0.0, 0.5, string1);
set(h,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle',...
    'Interpreter','none',...
    'FontSize',8)

switch lower(main.trusted)
    case 'true'
        h = text(0.0, 1.0,'TRUSTED OBSERVATION. Values include outlying obs. ');
        set(h,'HorizontalAlignment','center', ...
            'VerticalAlignment','middle',...
            'Interpreter','none',...
            'FontSize',20)
    otherwise
end

%=====================================================================


function [y,ydims] = FindTemporalVars(x)
%% Returns UNIQUE (i.e. base) temporal variable names

global verbose

if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
    error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;

for i = 1:length(x.allvarnames)
    indx = strfind(x.allvardims{i},'time');
    if (indx > 0)
        j = j + 1;
        basenames{j} = ReturnBase(x.allvarnames{i}); %#ok<AGROW>
        basedims{ j} = x.allvardims{i}; %#ok<AGROW>
    end
end

[~,i,~] = unique(basenames);
y     = cell(length(i),1);
ydims = cell(length(i),1);
for k = 1:length(i)
    if (verbose), fprintf('%3d is %s\n',k,basenames{i(k)}); end
    y{k}     = basenames{i(k)};
    ydims{k} = basedims{ i(k)};
end


%=====================================================================


function s = ReturnBase(string1)
%% Pick off the variable name.
inds = strfind(string1,'_guess');
if (inds > 0 )
    s = string1(1:inds-1);
    return
end
inds = strfind(string1,'_analy');
if (inds > 0 )
    s = string1(1:inds-1);
    return
end
inds = strfind(string1,'_VPguess');
if (inds > 0 )
    s = string1(1:inds-1);
    return
end
inds = strfind(string1,'_VPanaly');
if (inds > 0 )
    s = string1(1:inds-1);
    return
end


%=====================================================================


function x = FindRange(y)
%% Trying to pick 'nice' limits for plotting.
% Completely ad hoc ... and not well posed.
%
% In this scope, y is bounded from below by 0.0
%
% If the numbers are very small ...

bob  = [y.ges_copy(:) ; ...
    y.anl_copy(:)];
inds = find(isfinite(bob));

if ( isempty(inds) )
    x = [0 1];
else
    glommed = bob(inds);
    ymin    = min(glommed);
    ymax    = max(glommed);
    
    if ( ymax > 1.0 )
        ymin = floor(min(glommed));
        ymax =  ceil(max(glommed));
    elseif ( ymax < 0.0 && strcmp(y.copystring,'bias') )
        ymax = 0.0;
    end
    
    Yrange = [ymin ymax];
    
    x = [min([Yrange(1) 0.0]) Yrange(2)];
end


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

%=====================================================================

function [h, legstr] = plot_quantity(phase, plotdat)

global figuredata

switch lower(phase)
    case 'prior'
        data      = squeeze(plotdat.ges_copy( plotdat.region,:,:,:));
        Nused     = squeeze(plotdat.ges_Nused(plotdat.region,:,:,:));
        color     = figuredata.ges_color;
        marker    = figuredata.marker1;
        linestyle = figuredata.solid;
        linewidth = figuredata.linewidth;
        string1   = 'forecast:';
    case 'posterior'
        data      = squeeze(plotdat.anl_copy( plotdat.region,:,:,:));
        Nused     = squeeze(plotdat.anl_Nused(plotdat.region,:,:,:));
        color     = figuredata.anl_color;
        marker    = figuredata.marker2;
        linestyle = figuredata.solid;
        linewidth = figuredata.linewidth;
        string1   = 'analysis:';
    otherwise
        error('phase (%s) not supported',phase)
end

% Determine legend text
if sum(Nused(:)) > 1
    data_mean = mean(data(isfinite(data)));
    legstr = sprintf('%s mean = %.5g', string1, data_mean);
else
    legstr = ' ';
end

h = line(plotdat.bincenters,data);
set(h, 'LineStyle',    linestyle, ...
    'LineWidth',       linewidth, ...
    'Color',           color, ...
    'Marker',          marker, ...
    'MarkerFaceColor', color, ...
    'MarkerSize', figuredata.MarkerSize);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
