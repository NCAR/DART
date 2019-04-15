function plotdat = plot_profile(fname, copy, varargin)
%% plot_profile plots the vertical profile of the observation-space quantities for all possible levels, all possible variables.
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
% USAGE: plotdat = plot_profile(fname, copy [,varargin]);
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
%
% range    : 'range' of the value being plotted. Default is to
%                automatically determine range based on the data values.
%
% verbose  : true/false to control amount of run-time output
%
% MarkerSize  : integer controlling the size of the symbols
%
% pause  : true/false to conrol pausing after each figure is created.
%          true will require hitting any key to continue to next plot
%
% OUTPUT: 'plotdat' is a structure containing what was plotted.
%         A .pdf of each graphic is created. Each .pdf has a name that
%         reflects the variable, quantity, and region being plotted.
%
% EXAMPLE 1: All the observation types possible are plotted in separate figures.
%
% fname   = 'obs_diag_output.nc';
% copy    = 'totalspread';
% plotdat = plot_profile(fname, copy);
%
% EXAMPLE 2: Just a single observation type.
%
% fname   = 'obs_diag_output.nc';
% copy    = 'totalspread';
% obsname = 'RADIOSONDE_TEMPERATURE';
% plotdat = plot_profile(fname, copy, 'obsname', obsname);
%
% EXAMPLE 3: You can tell it to be quiet, too.
% plotdat = plot_profile(fname, copy, 'verbose', 'no');

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------

default_obsname    = 'none';
default_verbosity  = true;
default_markersize = 12;
default_pause      = false;
default_range      = [NaN NaN];
p = inputParser;

addRequired(p,'fname',@ischar);
addRequired(p,'copy',@ischar);
if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'obsname',    default_obsname,    @ischar);
    addParameter(p,'verbose',    default_verbosity,  @islogical);
    addParameter(p,'MarkerSize', default_markersize, @isnumeric);
    addParameter(p,'pause',      default_pause,      @islogical);
    addParameter(p,'range',      default_range,      @isnumeric);
else
    addParamValue(p,'obsname',   default_obsname,    @ischar);    %#ok<NVREPL>
    addParamValue(p,'verbose',   default_verbosity,  @islogical); %#ok<NVREPL>
    addParamValue(p,'MarkerSize',default_markersize, @isnumeric); %#ok<NVREPL>
    addParamValue(p,'pause',     default_pause,      @islogical); %#ok<NVREPL>
    addParamValue(p,'range',     default_range,      @isnumeric); %#ok<NVREPL>
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

% Either use all the variables or just the one optionally specified.
if (nvars == 0)
    [plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(fname);
    [plotdat.varnames,    plotdat.vardims]    = FindVerticalVars(plotdat);
    plotdat.nvars       = length(plotdat.varnames);
else
    plotdat.varnames{1} = obsname;
    plotdat.nvars       = nvars;
end

global figuredata verbose

figuredata            = set_obsdiag_figure('tall');
figuredata.MarkerSize = p.Results.MarkerSize;
verbose               = p.Results.verbose;

%%---------------------------------------------------------------------
% Loop around (copy-level-region) observation types
%----------------------------------------------------------------------

for ivar = 1:plotdat.nvars
    
    % create the variable names of interest.
    
    plotdat.myvarname = plotdat.varnames{ivar};
    plotdat.guessvar  = sprintf('%s_VPguess',plotdat.varnames{ivar});
    plotdat.analyvar  = sprintf('%s_VPanaly',plotdat.varnames{ivar});
    
    plotdat.trusted   = nc_read_att(fname, plotdat.guessvar, 'TRUSTED');
    if (isempty(plotdat.trusted)), plotdat.trusted = 'NO'; end
    
    % get appropriate vertical coordinate variable
    
    [dimnames,~] = nc_var_dims(fname, plotdat.guessvar);
    
    % this is a superfluous check ... FindVerticalVars already weeds out
    % variables only present on surface or undef because obs_diag
    % does not time-average statistics for these.
    
    if (~ isempty(strfind(dimnames{2},'surface')))
        fprintf('%s is a surface field.\n',plotdat.guessvar)
        fprintf('Cannot display a surface field this way.\n')
        continue
    elseif (~ isempty(strfind(dimnames{2},'undef')))
        fprintf('%s has no vertical definition.\n',plotdat.guessvar)
        fprintf('Cannot display this field this way.\n')
        continue
    end
    
    [levels, level_units, nlevels, level_edges, Yrange] = FindVerticalInfo(fname, plotdat.guessvar);
    plotdat.levels      = levels;
    plotdat.level_units = level_units;
    plotdat.nlevels     = nlevels;
    plotdat.level_edges = level_edges;
    plotdat.Yrange      = Yrange;
    
    % Matlab likes strictly ASCENDING order for the axes and ticks,
    % then you can impose the direction.
    
    if (plotdat.levels(1) > plotdat.levels(plotdat.nlevels))
        plotdat.YDir = 'reverse';
    else
        plotdat.YDir = 'normal';
    end
    
    [levels, ~]   = sort(plotdat.levels);
    plotdat.YTick = unique(levels);
    
    % Add error-checking for output from older versions of obs_diag.
    if (length(plotdat.YTick) ~= length(plotdat.levels))
        error('There is a duplicated value in the array specifying the levels - must change your input.nml and rerun obs_diag')
    end
    
    level_edges         = sort(plotdat.level_edges);
    plotdat.level_edges = level_edges;
    
    % guess(nregions,nlevels,ncopies)
    
    guess = ncread(fname, plotdat.guessvar);
    analy = local_ncread(fname, plotdat.analyvar);
    if (isempty(analy))
        analy    = guess;  % make the variable the same shape as guess
        analy(:) = NaN;    % and fill it with nothing
        plotdat.has_analysis = false;
        plotdat.post_string = '';
    else
        plotdat.has_analysis = true;
        plotdat.post_string = '; \diamondsuit=posteriorOK';
    end
    
    % check to see if there is anything to plot
    % The number possible is decreased by the number of observations
    % rejected by namelist control.
    
    priorQCs = get_qc_values(fname, plotdat.guessvar, ...
        'fatal', false, ...
        'verbose', verbose);
    
    plotdat.ges_Neval = priorQCs.num_evaluated;
    plotdat.ges_Nposs = priorQCs.nposs;
    plotdat.ges_Nused = priorQCs.nused;
    plotdat.ges_copy  = guess(:,:,plotdat.copyindex);
    
    if ( sum(plotdat.ges_Nposs(:)) < 1 )
        fprintf('no obs for %s...  skipping\n', plotdat.varnames{ivar})
        continue
    end
    
    if (plotdat.has_analysis)
        posteQCs = get_qc_values(fname, plotdat.analyvar, ...
            'fatal', false, ...
            'verbose', verbose);
        plotdat.anl_Nused = posteQCs.nused;
        plotdat.anl_copy  = analy(:,:,plotdat.copyindex);
    else
        plotdat.anl_Nused = zeros(size(plotdat.ges_Nused));
        plotdat.anl_copy  = plotdat.ges_copy;  % needed for determining limits
    end
    
    % call report_qc_values.m
    
    plotdat.Xrange     = FindRange(plotdat);
    
    % plot by region - each in its own figure.
    
    for iregion = 1:plotdat.nregions
        figure(iregion);
        clf(iregion);
        orient(figuredata.orientation);
        plotdat.region   = iregion;
        plotdat.myregion = deblank(plotdat.region_names(iregion,:));
        
        myplot(plotdat);
        
        BottomAnnotation(fname)
        
        psfname = sprintf('%s_%s_profile_region%d', ...
            plotdat.varnames{ivar}, plotdat.copystring, iregion);
        
        if verLessThan('matlab','R2016a')
            print(gcf, '-dpdf', psfname);
        else
            print(gcf, '-dpdf', '-bestfit', psfname);
        end
        
        % block to go slow and look at each one ...
        if (p.Results.pause)
            disp('Pausing, hit any key to continue ...')
            pause
        end
        
    end
end

%=====================================================================
% 'Helper' functions
%=====================================================================


function myplot(plotdat)

%% The prior and posterior are plotted as separate items.

% Plot 'xxx' on the bottom axis.
% The observation count will use the axis on the top.
% Ultimately, we want to suppress the 'auto' feature of the
% axis labelling, so we manually set some values that normally
% don't need to be set.

global figuredata verbose

ax1 = subplot('position',figuredata.position);
orient(figuredata.orientation)

% add type of vertical coordinate info for adjusting axes to accomodate legend

Stripes(plotdat.Xrange, plotdat.level_edges, plotdat.level_units);
set(ax1, 'YDir', plotdat.YDir, 'YTick', plotdat.YTick, 'Layer', 'top')
set(ax1,'YAxisLocation','left','FontSize',figuredata.fontsize)

% draw the result of the experiment

[hprior, prior_legstr] = plot_quantity('prior', plotdat);

ges_Nposs = plotdat.ges_Nposs(plotdat.region,:);
ges_Nused = plotdat.ges_Nused(plotdat.region,:);
anl_Nused = plotdat.anl_Nused(plotdat.region,:);
anl_Ngood = sum(anl_Nused);

if verbose
    fprintf('region %d %s nobs_poss %d prior %d poste %d\n', ...
        plotdat.region, plotdat.myvarname, ...
        sum(ges_Nposs), sum(ges_Nused), anl_Ngood)
    fprintf('region %d %s\n\n', plotdat.region, prior_legstr)
end

if anl_Ngood > 0
    [hposte, poste_legstr] = plot_quantity('posterior', plotdat);
    h = legend([hprior, hposte], prior_legstr, poste_legstr);
else
    h = legend(hprior,prior_legstr);
end

set(h,'Interpreter','none','Box','off','Location','NorthWest')

if verLessThan('matlab','R2017a')
    % Convince Matlab to not autoupdate the legend with each new line.
    % Before 2017a, this was the default behavior, so do nothing.
    % We do not want to add the bias line to the legend, for example.
else
    h.AutoUpdate = 'off';
end

switch plotdat.copystring
    case {'bias'}
        zeroline = line([0 0],plotdat.Yrange,'Color',[200 200 200]/255,'Parent',ax1);
        set(zeroline,'LineWidth',2.5,'LineStyle','-')
        plotdat.xlabel = sprintf('%s (%s)',plotdat.copystring, plotdat.biasconv);
    otherwise
        plotdat.xlabel = plotdat.copystring;
end

% If the observation is trusted, reference that somehow

switch lower(plotdat.trusted)
    case 'true'
        axlims = axis;
        tx = axlims(2) + (axlims(2) - axlims(1))/20;
        if  strcmpi('normal',plotdat.YDir)
            ty = plotdat.Yrange(1);
        else
            ty = plotdat.Yrange(2);
        end
        h = text(tx,ty,'TRUSTED. Values include outlying observations.');
        set(h,'FontSize',20,'Rotation',90,'VerticalAlignment','middle')
    otherwise
end

% Create another axes to use for plotting the observation counts

ax2 = axes('position',get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'YAxisLocation','right', ...
    'Color','none', ...
    'XColor',figuredata.obs_color, ...
    'YColor',get(ax1,'YColor'), ...
    'YLim',get(ax1,'YLim'), ...
    'YDir',get(ax1,'YDir'), ...
    'FontSize',get(ax1,'FontSize'));

ax2h1 = line(ges_Nposs, plotdat.levels, 'Parent', ax2);
ax2h2 = line(ges_Nused, plotdat.levels, 'Parent', ax2);

set(ax2h1, 'LineStyle', 'none', ...
    'Color',      figuredata.obs_color, ...
    'Marker',     figuredata.obs_marker, ...
    'MarkerSize', figuredata.MarkerSize);

set(ax2h2, 'LineStyle', 'none', ...
    'Color',      figuredata.ges_color, ...
    'Marker',     figuredata.ges_marker, ...
    'MarkerSize', figuredata.MarkerSize);

if anl_Ngood > 0
    ax2h3 = line(anl_Nused, plotdat.levels, 'Parent',ax2);
    set(ax2h3, 'LineStyle', 'none', ...
        'Color',     figuredata.anl_color, ...
        'Marker',    figuredata.anl_marker, ...
        'MarkerSize',figuredata.MarkerSize);
end

% use same Y ticks - but no labels.
set(ax2,'YTick',get(ax1,'YTick'), 'YTicklabel',[]);

% use the same X ticks, but find the right label values
xscale = matchingXticks(ax1,ax2);

set(get(ax1,'Ylabel'),'String',plotdat.level_units, ...
    'Interpreter','none','FontSize',figuredata.fontsize)
set(get(ax1,'Xlabel'),'String',{plotdat.xlabel, plotdat.timespan}, ...
    'Interpreter','none','FontSize',figuredata.fontsize)

% determine if the observation was flagged as 'evaluate' or 'assimilate'

if sum(plotdat.ges_Neval(plotdat.region,:)) > 0
    string1 = sprintf('# of obs (o=possible; %s %s) x %d', ...
        '\ast=evaluated', plotdat.post_string, uint32(xscale));
else
    string1 = sprintf('# of obs (o=possible; %s %s) x %d', ...
        '\ast=assimilated', plotdat.post_string, uint32(xscale));
end

set(get(ax2,'Xlabel'), 'String', string1, 'FontSize', figuredata.fontsize)

title({plotdat.myregion, plotdat.myvarname},  ...
    'Interpreter', 'none', 'FontSize', figuredata.fontsize, 'FontWeight', 'bold')


%=====================================================================


function BottomAnnotation(main)
%% annotates the filename containing the data being plotted
subplot('position',[0.10 0.01 0.8 0.04])
axis off

if ( main(1) == '/' )   % must be an absolute pathname
    string1 = sprintf('data file: %s',main);
else
    mydir = pwd;
    string1 = sprintf('data file: %s/%s',mydir,main);
end

h = text(0.5, 0.33, string1);
set(h, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment','middle', ...
    'Interpreter', 'none', ...
    'FontSize', 8);


%=====================================================================


function [y,ydims] = FindVerticalVars(x)
%% Returns UNIQUE (i.e. base) vertical variable names
% In this context, if the variable has a 'time' dimension
% it cannot be a variable of interest.

global verbose

if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
    error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;
basenames = struct([]);
basedims  = struct([]);

for i = 1:length(x.allvarnames)
    dimnames = lower(x.allvardims{i});
    if (isempty(strfind(dimnames,'time'))) %#ok<STREMP>
        platform = ReturnBase(x.allvarnames{i});
        if (~ isempty(platform))
            j = j + 1;
            basenames{j} = platform;
            basedims{j}  = x.allvardims{i};
        end
    end
end

[~,i,~] = unique(basenames);
y       = struct([]);
ydims   = struct([]);

for k = 1:length(i)
    if (verbose), fprintf('%3d is %s\n',k,basenames{i(k)}); end
    y{k} = basenames{i(k)};
    ydims{k} = basedims{i(k)};
end


%=====================================================================


function [levels, level_units, nlevels, level_edges, Yrange] = FindVerticalInfo(fname,varname)
%% Find the vertical dimension and harvest some info

varinfo  = ncinfo(fname,varname);
leveldim = [];

for i = 1:length(varinfo.Dimensions)
    inds = strfind(varinfo.Dimensions(i).Name,'level');
    if ( ~ isempty(inds)), leveldim = i; end
end

if ( isempty(leveldim) )
    error('There is no level information for %s in %s',varname,fname)
end

levels      = ncread(   fname,varinfo.Dimensions(leveldim).Name);
level_units = ncreadatt(fname,varinfo.Dimensions(leveldim).Name,'units');
nlevels     = varinfo.Size(leveldim);
edgename    = sprintf('%s_edges',varinfo.Dimensions(leveldim).Name);
level_edges = ncread(fname, edgename);
Yrange      = [min(level_edges) max(level_edges)];


%=====================================================================


function s = ReturnBase(string1)
%% Pick off the variable name.
s = [];
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

bob  = [y.ges_copy(:); y.anl_copy(:)];
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
    end
    
    if (ymin == 0 && ymax == 0)
        ymax = 1;
    end
    
    if (ymin == ymax)
        ymin = ymin - 0.1*ymin;
        ymax = ymax + 0.1*ymax;
    end
    
    Yrange = [ymin ymax];
    
    % Make sure a zero bias is visible on plot
    switch y.copystring
        case {'bias'}
            if  ymax < 0
                Yrange = [ ymin 0.0 ];
            elseif  ymin > 0
                Yrange = [ 0.0 ymax ];
            end
    end
    
    x = sort([min([Yrange(1) 0.0]) Yrange(2)] ,'ascend');
end


%=====================================================================


function h = Stripes(x,edges,units)
%% plot the subtle background stripes that indicate the vertical
%  extent of what was averaged.
%
% FIXME:
% This really should be modified to add a percentage of the data
% range to provide space for the legend. Right now it is hardwired
% to assume that we are plotting hPa, on a 'reverse' axis.
% kdr axlims(3) should be conditional on the observation vertical coordinate:
%     values for pressure coordinates are inappropriate for height coord.
%     It also assumes 4 plots/page, but 2 works better for plotting all levels of CAM5.
%     That requires a smaller % of vertical range for the legend.

% plot two little dots at the corners to make Matlab
% choose some plot limits. Given those nice limits and
% tick labels ... KEEP THEM. Later, make the dots invisible.

h = plot([min(x) max(x)],[min(edges) max(edges)]);
axlims          = axis;
legend_fraction = 0.1375;

% partial fix to legend space; add in option for vert coord = height.

switch lower(units)
    case 'hpa'
        axlims(4) = max(edges);
        axlims(3) = min(edges) - legend_fraction*(axlims(4)-min(edges));
    case 'm'
        axlims(3) = min(edges);
        axlims(4) = max(edges) + legend_fraction*(max(edges)-axlims(3));
    otherwise
end
axis(axlims)

% set up list of x,y values defining corner of every other stripe.

xc = [ axlims(1) axlims(2) axlims(2) axlims(1) axlims(1) ];

hold on;
for i = 1:2:(length(edges)-1)
    yc = [ edges(i) edges(i) edges(i+1) edges(i+1) edges(i) ];
    fill(xc,yc,[0.8 0.8 0.8],'EdgeColor','none');
end
hold off;

set(gca,'XGrid','on')
set(h,'Visible','off') % make the dots invisible


%=====================================================================


function value = local_ncread(fname,varname)
%% If the variable exists in the file, return the contents of the variable.
% if the variable does not exist, return empty value instead of error-ing
% out.

[variable_present, ~] = nc_var_exists(fname,varname);
if (variable_present)
    value = ncread(fname, varname);
else
    value = [];
end

%=====================================================================

function [h, legstr] = plot_quantity(phase, plotdat)

global figuredata

switch lower(phase)
    case 'prior'
        data      = plotdat.ges_copy( plotdat.region,:);
        Nused     = plotdat.ges_Nused(plotdat.region,:);
        color     = figuredata.ges_color;
        marker    = figuredata.marker1;
        linestyle = figuredata.solid;
        linewidth = figuredata.linewidth;
        string1   = 'forecast:';
    case 'posterior'
        data      = plotdat.anl_copy( plotdat.region,:);
        Nused     = plotdat.anl_Nused(plotdat.region,:);
        color     = figuredata.anl_color;
        marker    = figuredata.marker2;
        linestyle = figuredata.solid;
        linewidth = figuredata.linewidth;
        string1   = 'analysis:';
    otherwise
        error('phase (%s) not supported',phase)
end

% Determine legend text
nobs = sum(Nused);
if ( nobs > 1 )
    data_mean = mean(data(isfinite(data)));
    legstr = sprintf('%s mean= %.5g', string1, data_mean);
else
    legstr = ' ';
end

h = line(data, plotdat.levels);
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
