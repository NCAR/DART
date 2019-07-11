function two_experiments_profile(files, titles, obsnames, copy, prpo, varargin)
% Plot two or more experiments on the same axis.
% Part of the observation-space diagnostics routines.
%
% 'obs_diag' produces a netcdf file containing the diagnostics.
% obs_diag condenses the obs_seq.final information into summaries for a few specified
% regions - on a level-by-level basis.
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
% USAGE: two_experiments_profile(files, titles, obsnames, copy, prpo, 'level', 1)
%
% files    : Cell array containing the locations of the obs_diag_output.nc
%            files to compare. Each file is presumed to be the results from
%            a single experiment.
%
% titles   : Cell array containing the titles used to annotate each of the experiments.
%
% obsnames : Cell array containing the strings of each observation type to plot.
%            Each observation type will be plotted in a separate graphic.
%
% copy     : string defining the metric of interest. 'rmse', 'spread', etc.
%            Possible values are available in the netcdf 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% prpo     : string defining whether to plot the prior or posterior metrics.
%            Due to the amount of information already plotted, we made a
%            conscious decision not to support plotting both prior and posterior
%            on the same plot.
%
% varargin: optional parameter-value pairs. Supported parameters are described below.
%
% level    : The index of the level to plot. Defaults to level 1.
%
% verbose  : true/false to control amount of run-time output
%
% MarkerSize  : integer controlling the size of the symbols
%
% pause    : true/false to conrol pausing after each figure is created.
%            true will require hitting any key to continue to next plot
%
% range    : 'range' of the value being plotted. Default is to
%            automatically determine range based on the data values.
%
% OUTPUT: A .pdf of each graphic is created. Each .pdf has a name that
%         reflects the variable, quantity, and region being plotted.
%
% EXAMPLE
%
% files = {'/ptmp/nancy/CSL/Base5/032-061s0_def_reg/obs_diag_output.nc',
%          '/ptmp/thoar/GPS+AIRS/Sep_032-061/obs_diag_output.nc'};
% titles = {'Base5', 'GPS+AIRS'};
% files = {'/glade/scratch/nancy/fvdiags/obs_diag_2005_08.nc', ...
%  '/glade/scratch/raeder/SE_NCEP_assim1/Diag_hemi_poles_2005.8.1-30/obs_diag_output.nc'};
% titles   = {'FV', 'SE'};
% obsnames = {'RADIOSONDE_TEMPERATURE', ...
%             'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_V_WIND_COMPONENT'};
% copy     = 'rmse';     % rmse, spread, totalspread, bias, etc.
% prpo     = 'analysis'; % [analy, analysis, posterior ] == posterior
% prpo     = 'forecast'; % [guess, forecast, prior     ] == prior
% prpo     = 'both';
%
% two_experiments_profile(files, titles, obsnames, copy, prpo)
%
% Example 2: restrict the data limits to the data in a certain vertical area.
%            In this case, the observations using a pressure vertical coordinate
%            between +Inf (the surface) and 100hPa (inclusive) are used to determine
%            the scale. All values will be plotted, but the highest levels may be
%            clipped. The optional argument pairs at the end consist of a string
%            and a length 2 array specifying the [bottom top] levels to consider.
%
% two_experiments_profile(files, titles, obsnames, copy, prpo,'plevel',[Inf 100])
%
% two_experiments_profile(files, titles, obsnames, copy, prpo, ...
%            'plevel',[Inf 100],'mlevel',[1 10],'hlevel',[0 20000])

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------
defaultPlevels = [ Inf  0 ];
defaultHlevels = [-Inf Inf];
defaultMlevels = [  1  Inf];

default_verbosity  = true;
default_markersize = 12;
default_pause      = false;
default_range      = [NaN NaN];
default_level      = 1;
p = inputParser;
addRequired(p,'files',@iscell);
addRequired(p,'titles',@iscell);
addRequired(p,'obsnames',@iscell);
addRequired(p,'copy',@ischar);
addRequired(p,'prpo',@ischar);

if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'plevel',defaultPlevels,@isnumeric);
    addParameter(p,'hlevel',defaultHlevels,@isnumeric);
    addParameter(p,'mlevel',defaultMlevels,@isnumeric);
    addParameter(p,'verbose',    default_verbosity,  @islogical);
    addParameter(p,'MarkerSize', default_markersize, @isnumeric);
    addParameter(p,'pause',      default_pause,      @islogical);
    addParameter(p,'range',      default_range,      @isnumeric);
    addParameter(p,'level',      default_level,      @isnumeric);
else
    addParamValue(p,'plevel',defaultPlevels,@isnumeric); %#ok<NVREPL>
    addParamValue(p,'hlevel',defaultHlevels,@isnumeric); %#ok<NVREPL>
    addParamValue(p,'mlevel',defaultMlevels,@isnumeric); %#ok<NVREPL>
    addParamValue(p,'verbose',   default_verbosity,  @islogical); %#ok<NVREPL>
    addParamValue(p,'MarkerSize',default_markersize, @isnumeric); %#ok<NVREPL>
    addParamValue(p,'pause',     default_pause,      @islogical); %#ok<NVREPL>
    addParamValue(p,'range',     default_range,      @isnumeric); %#ok<NVREPL>
    addParamValue(p,'level',     default_level,      @isnumeric); %#ok<NVREPL>
end

p.parse(files, titles, obsnames, copy, prpo, varargin{:});

% if you want to echo the input
% disp(['files   : ', p.Results.files])
% disp(['titles  : ', p.Results.titles])
% disp(['obsnames: ', p.Results.obsnames])
% disp(['copy    : ', p.Results.copy])
% disp(['prpo    : ', p.Results.prpo])
% fprintf( 'plevel : %f %f \n', p.Results.plevel)
% fprintf( 'hlevel : %f %f \n', p.Results.hlevel)
% fprintf( 'mlevel : %f %f \n', p.Results.mlevel)

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

% if ~isempty(p.UsingDefaults)
%    disp('Using defaults: ')
%    disp(p.UsingDefaults)
% end

if (numel(p.Results.plevel) ~= 2)
    error('plevel must be an array of length two ... [bottom top]')
end

if (numel(p.Results.hlevel) ~= 2)
    error('hlevel must be an array of length two ... [bottom top]')
end

if (numel(p.Results.mlevel) ~= 2)
    error('mlevel must be an array of length two ... [bottom top]')
end

% Now that the input passes sanity checks ...

if (length(files) ~= length(titles))
    error('each file must have a title')
end

NumExp = length(files);

for i = 1:NumExp
    if (exist(files{i},'file') ~= 2)
        error('File %s does not exist',files{i})
    end
end

if (NumExp ~= length(titles))
    error('each file must have an experiment title')
end

%% set up all the stuff that is common.

global figuredata

commondata = check_compatibility(files, prpo, obsnames, copy);
figuredata = set_obsdiag_figure('tall', 'numexp', NumExp);
figuredata.MarkerSize = p.Results.MarkerSize;

%%--------------------------------------------------------------------
% Set some static data
%---------------------------------------------------------------------

nvars = length(obsnames);

for ivar = 1:nvars
    fprintf('Working on %s ...\n',obsnames{ivar})
    
    for iregion = 1:commondata.nregions
        
        %---------------------------------------------------------------------
        % Getting the data for each experiment
        %---------------------------------------------------------------------
        
        Dlimits = zeros(NumExp,2);  % range of the data
        plotobj = cell(1,NumExp);
        
        for iexp = 1:NumExp
            
            plotobj{iexp} = getvals(files{iexp}, commondata.targets{ivar}, copy, iregion, p);
            plotobj{iexp}.title  = titles{iexp};
            plotobj{iexp}.phase  = commondata.phase;
            
            Dlimits(iexp,:) = plotobj{iexp}.Drange;
            
        end
        
        %---------------------------------------------------------------------
        % Find nice limits that encompass all experiments
        % Note that Dlimits has been constructed by ignoring the top levels.
        %---------------------------------------------------------------------
        
        Drange = [min(Dlimits(:,1)) max(Dlimits(:,2))];
        span = abs(Drange(2) - Drange(1))* 0.05;
        Drange(1) = Drange(1) - span;
        Drange(2) = Drange(2) + span;
        
        %---------------------------------------------------------------------
        % Plot all regions - one region to a page
        %---------------------------------------------------------------------
        
        myplot(plotobj, Drange);
        
        BottomAnnotation(files)
        
        psfname = sprintf('%s_%s_region%d_profile_%dexp', ...
            obsnames{ivar}, plotobj{1}.copystring, iregion, NumExp);
        
        if verLessThan('matlab','R2016a')
            print(iregion, '-dpdf', psfname)
        else
            print(iregion, '-dpdf', '-bestfit', psfname)
        end
        
    end % of loop around regions
    
    if ( ivar ~= nvars && p.Results.pause )
        disp('Pausing, hit any key to continue ...')
        pause
    end
    
end  % of loop around variable



%=====================================================================
% End of main function body. Helper functions below.
%=====================================================================



function common = check_compatibility(filenames, prpo, varnames, copystring)
%% Trying to prevent the comparison of apples and oranges.
% make sure the diagnostics were generated the same way.

% need to check that the timeframe is the same for all files

% need to check that the region definitions are the same for all files

mystat     = 0;
nexp       = length(filenames);
commondata = cell(1,nexp);
targets    = struct([]);

for i = 1:length(varnames)
    switch lower(prpo)
        case {'guess','forecast','prior'}
            targets{i} = sprintf('%s_VPguess',varnames{i});
            commondata{i}.phase = 'prior';
        case {'analy','analysis','posterior'}
            targets{i} = sprintf('%s_VPanaly',varnames{i});
            commondata{i}.phase = 'posterior';
        otherwise
            error('unknown prpo ... "%s"',prpo)
    end
end

for i = 1:nexp
    
    varexist(filenames{i}, {targets{:}, 'time', 'time_bounds'}) %#ok<CCAT>
    
    commondata{i}.targets      = targets;
    commondata{i}.region_names = strtrim(ncread(filenames{i},'region_names')');
    commondata{i}.times        = ncread(filenames{i}, 'time');
    commondata{i}.time_bnds    = ncread(filenames{i}, 'time_bounds');
    commondata{i}.copyindex    = get_copy_index(filenames{i},copystring);
    commondata{i}.nobstypes    = nc_dim_info(filenames{i}, 'obstypes');
    commondata{i}.nregions     = nc_dim_info(filenames{i}, 'region');
    commondata{i}.time_to_skip = nc_read_att(filenames{i}, '/', 'time_to_skip');
    commondata{i}.lonlim1      = nc_read_att(filenames{i}, '/', 'lonlim1');
    commondata{i}.lonlim2      = nc_read_att(filenames{i}, '/', 'lonlim2');
    commondata{i}.latlim1      = nc_read_att(filenames{i}, '/', 'latlim1');
    commondata{i}.latlim2      = nc_read_att(filenames{i}, '/', 'latlim2');
end

% error checking - compare everything to the first experiment
for i = 2:nexp
    
    if (any(commondata{i}.lonlim1 ~= commondata{1}.lonlim1))
        fprintf('The left longitudes of the regions (i.e. lonlim1) are not compatible.\n')
        mystat = 1;
    end
    
    if (any(commondata{i}.lonlim2 ~= commondata{1}.lonlim2))
        fprintf('The right longitudes of the regions (i.e. lonlim2) are not compatible.\n')
        mystat = 1;
    end
    
    if (any(commondata{i}.latlim1 ~= commondata{1}.latlim1))
        fprintf('The bottom latitudes of the regions (i.e. latlim1) are not compatible.\n')
        mystat = 1;
    end
    
    if (any(commondata{i}.latlim2 ~= commondata{1}.latlim2))
        fprintf('The top latitudes of the regions (i.e. latlim2) are not compatible.\n')
        mystat = 1;
    end
    
    if (any(commondata{i}.time_bnds ~= commondata{1}.time_bnds))
        fprintf('The time boundaries of the experiments (i.e. time_bnds) are not compatible.\n')
        mystat = 1;
    end
    
    if (any(commondata{i}.time_to_skip ~= commondata{1}.time_to_skip))
        fprintf('The time skipped in the experiments (i.e. time_to_skip) are not compatible.\n')
        mystat = 1;
    end
    
end

if mystat > 0
    error('The experiments are not compatible ... stopping.')
end

common = commondata{1};

% Coordinate between time types and dates

timeunits         = nc_read_att(filenames{1},'time','units');
timebase          = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin        = datenum(timebase(1),timebase(2),timebase(3));
timefloats        = zeros(size(commondata{1}.time_to_skip));  % stupid int32 type conversion
timefloats(:)     = commondata{1}.time_to_skip(:);
skip_seconds      = timefloats(4)*3600 + timefloats(5)*60 + timefloats(6);
iskip             = timefloats(3) + skip_seconds/86400;

common.bincenters = commondata{1}.times     + timeorigin;
common.binedges   = commondata{1}.time_bnds + timeorigin;
common.Nbins      = length(common.bincenters);
common.toff       = common.binedges(1) + iskip;

common.timespan   = sprintf('%s through %s', datestr(common.toff), ...
    datestr(max(common.binedges(:))));


%=====================================================================


function plotdat = getvals(fname, varname, copystring, regionindex, opt )
%% basic function to retrieve plotting data

if (exist(fname,'file') ~= 2)
    error('%s does not exist',fname)
end

plotdat            = read_obsdiag_staticdata(fname,copystring);
plotdat.region     = regionindex;
plotdat.varname    = varname;

myinfo.diagn_file  = fname;
myinfo.copyindex   = plotdat.copyindex;
myinfo.regionindex = plotdat.region;
[start, count]     = GetNCindices(myinfo,'diagn',plotdat.varname);
hyperslab          = ncread(fname, plotdat.varname, start, count);
plotdat.data       = squeeze(hyperslab);
plotdat.trusted    = nc_read_att(fname, plotdat.varname, 'TRUSTED');
if (isempty(plotdat.trusted)), plotdat.trusted = 'NO'; end

% Now that we know the variable ... get the appropriate vertical information

varinfo               = ncinfo(fname,plotdat.varname);
plotdat.levels        = ncread(fname,varinfo.Dimensions(2).Name);
plotdat.level_units   = nc_read_att(fname,varinfo.Dimensions(2).Name,'units');
plotdat.nlevels       = length(plotdat.levels);
plotdat.level_edges   = ncread(fname,sprintf('%s_edges',varinfo.Dimensions(2).Name));

plotdat.YDir = 'normal';
inds = 1:plotdat.nlevels;

%% Determine data limits
% find the levels of interest for setting the data limits

switch lower(varinfo.Dimensions(2).Name)
    case {'plevel'}
        plotdat.YDir = 'reverse';
        inds = find((plotdat.levels <= opt.Results.plevel(1)) & ...
            (plotdat.levels >= opt.Results.plevel(2)));
    case {'hlevel'}
        inds = find((plotdat.levels >= opt.Results.hlevel(1)) & ...
            (plotdat.levels <= opt.Results.hlevel(2)));
    case {'mlevel'}
        inds = find((plotdat.levels >= opt.Results.mlevel(1)) & ...
            (plotdat.levels <= opt.Results.mlevel(2)));
    otherwise
end

bob = plotdat.data(inds);

switch copystring
    case {'bias'}
        %  always make sure we have a zero bias line ...
        dmin = min( [ min(bob) 0.0 ] );
        dmax = max( [ max(bob) 0.0 ] );
        plotdat.Drange = [ dmin dmax ];
        plotdat.xlabel = sprintf('%s (%s)',copystring, plotdat.biasconv);
    case {'rmse'}
        plotdat.Drange = [0.0 max(bob)];
        plotdat.xlabel = copystring;
    otherwise
        plotdat.Drange = [min(bob) max(bob)];
        plotdat.xlabel = copystring;
end


qcvalues = get_qc_values(fname, plotdat.varname, ...
    'regionindex', plotdat.region, ...
    'fatal', false, ...
    'verbose', false);

plotdat.nposs         = qcvalues.nposs;
plotdat.nused         = qcvalues.nused;
plotdat.num_evaluated = qcvalues.num_evaluated;

if sum(plotdat.num_evaluated(:) > 0)
    plotdat.assim_eval_string = 'evaluated';
else
    plotdat.assim_eval_string = 'assimilated';
end

%=====================================================================


function myplot( plotdat, Drange)
%% Create graphic for one region - for all experiments.

global figuredata

Nexp    = length(plotdat);
iregion = plotdat{1}.region;

figure(iregion);
clf(iregion); orient(figuredata.orientation);
ax1 = subplot('position',figuredata.position);

Stripes(Drange, plotdat{1}.level_edges, plotdat{1}.level_units, Nexp);
set(ax1,'YDir',plotdat{1}.YDir,'YTick',sort(plotdat{1}.levels),'Layer','top')
set(ax1,'YAxisLocation','left','FontSize',figuredata.fontsize)

% draw the results of the experiments - each with their own line type.
hd     = [];   % handle to an unknown number of data lines
legstr = {[]}; % strings for the legend

for i = 1:Nexp
    hd(i) = line(plotdat{i}.data, plotdat{i}.levels, ...
        'Color',           figuredata.expcolors{i}, ...
        'Marker',          figuredata.expsymbols{i}, ...
        'MarkerSize',      figuredata.MarkerSize, ...
        'MarkerFaceColor', figuredata.expcolors{i}, ...
        'LineStyle',       figuredata.prpolines{1}, ...
        'LineWidth',       figuredata.linewidth,'Parent',ax1); %#ok<AGROW>
    
    legstr{i} = sprintf('%s %s',plotdat{i}.title,plotdat{i}.phase);
    
end

switch plotdat{1}.copystring
    case {'bias','rmse'}
        zeroline = line([0 0],get(ax1,'YLim'),'Color',[200 200 200]/255,'Parent',ax1);
        set(zeroline,'LineWidth',2.5,'LineStyle','-')
    otherwise
end

set(ax1,'XLim',Drange)

axlims = axis;
dx = (axlims(2) - axlims(1))/20;
if strcmpi('normal',get(ax1,'YDir'))
    ty = axlims(3);
else
    ty = axlims(4);
end

for i = 1:Nexp
    % If the observation is trusted, reference that somehow
    switch lower(plotdat{i}.trusted)
        case 'true'
            tx = axlims(2) + (i*dx);
            h = text(tx,ty,sprintf('TRUSTED OBSERVATION in %s',plotdat{i}.title));
            set(h, 'FontSize', 20, 'Rotation', 90, ...
                'VerticalAlignment', 'middle', 'Interpreter', 'none')
        otherwise
    end
end

% Create another axes to use for plotting the observation counts

ax2 = axes('position',get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'YAxisLocation','right',...
    'Color','none', ...
    'XColor','b', ...
    'YColor',get(ax1,'YColor'), ...
    'YLim',get(ax1,'YLim'), ...
    'YDir',get(ax1,'YDir'), ...
    'FontSize',get(ax1,'FontSize'));

% Plot the data, which sets the range of the axis
for i = 1:Nexp
    ax2h1 = line(plotdat{i}.nposs, plotdat{i}.levels, 'Parent', ax2);
    ax2h2 = line(plotdat{i}.nused, plotdat{i}.levels, 'Parent', ax2);
    
    set(ax2h1,'LineStyle','none', ...
        'Color',     figuredata.expcolors{i}, ...
        'Marker',    figuredata.obs_marker, ...
        'MarkerSize',figuredata.MarkerSize);
    
    set(ax2h2,'LineStyle','none', ...
        'Color',     figuredata.expcolors{i}, ...
        'Marker',    figuredata.ges_marker, ...
        'MarkerSize',figuredata.MarkerSize);
end

% use same Y ticks but no labels
set(ax2,'YTick',get(ax1,'YTick'), 'YTicklabel',[]);

% use the same X ticks, but find the right label values
xscale = matchingXticks(ax1,ax2);

% Annotate the whole thing - gets pretty complicated for multiple
% regions on one page. Trying to maximize content, minimize clutter.
% Any plot object will do for annotating region,levels,etc

annotate(ax1, ax2, plotdat{1}, xscale)

lh = legend(hd,legstr,'Location','NorthWest');
set(lh,'Interpreter','none','Box','off','FontSize',figuredata.fontsize);

if verLessThan('matlab','R2017a')
    % Convince Matlab to not autoupdate the legend with each new line.
    % Before 2017a, this was the default behavior, so do nothing.
    % We do not want to add the bias line to the legend, for example.
else
    lh.AutoUpdate = 'off';
end

%=====================================================================


function h = Stripes(x, edges, units, nexp)
%% plot the subtle background stripes that indicate the vertical
%  extent of what was averaged.
%
%  axlims(3) should be conditional on the observation vertical coordinate:
%  values for pressure coordinates are inappropriate for height coord.

% plot two little dots at the corners to make Matlab
% choose some plot limits. Given those nice limits and
% tick labels ... KEEP THEM. Later, make the dots invisible.

h = plot([min(x) max(x)],[min(edges) max(edges)]);
axlims          = axis;
legend_fraction = 0.05 * nexp + 0.005;

% partial fix to legend space; add in option for vert coord = height.

switch lower(units)
    case 'hpa'
        axlims(4) = max(edges);
        axlims(3) = min(edges) - legend_fraction*(axlims(4)-min(edges));
    case 'm'
        axlims(3) = min(edges) ;
        axlims(4) = max(edges) + legend_fraction*(max(edges)-axlims(3));
    otherwise
end
axis(axlims)

% set up list of x,y values defining corner of every other stripe

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


function annotate(ax1, ax2, plotobj, xscale)

%% One figure ... everything gets annotated.

global figuredata

set(get(ax1,'Ylabel'),'String',plotobj.level_units, ...
    'Interpreter','none','FontSize',figuredata.fontsize)
set(get(ax1,'Xlabel'),'String',{plotobj.xlabel,plotobj.timespan}, ...
    'Interpreter','none','FontSize',figuredata.fontsize)

string1 = sprintf('# of obs (o=possible, %s=%s ) x ', '\ast', ...
    plotobj.assim_eval_string);

set(get(ax2,'Xlabel'),'String', [string1 int2str(uint32(xscale))], ...
    'FontSize',figuredata.fontsize)

th = title({deblank(plotobj.region_names(plotobj.region,:)), plotobj.varname});
set(th,'Interpreter','none','FontSize',figuredata.fontsize,'FontWeight','bold');


%=====================================================================


function BottomAnnotation(filenames)
%% annotates the filenames containing the data being plotted

nfiles  = length(filenames);
dy      = 1.0/(nfiles+3);
yheight = 0.0225*(nfiles+1);

subplot('position',[0.10 0.01 0.8 yheight])
axis off

% list all the files.

for ifile = 1:nfiles
    main = filenames{ifile};
    
    fullname = which(main);   % Could be in MatlabPath
    if( isempty(fullname) )
        if ( main(1) == '/' )  % must be a absolute pathname
            string1 = sprintf('data file: %s',main);
        else                   % must be a relative pathname
            mydir = pwd;
            string1 = sprintf('data file: %s/%s',mydir,main);
        end
    else
        string1 = sprintf('data file: %s',fullname);
    end
    
    ty = 1.0 - (ifile+2)*dy;
    h = text(0.5, ty, string1);
    set(h, 'Interpreter', 'none', 'FontSize', 8);
    set(h, 'HorizontalAlignment','center');
end


%=====================================================================


function varexist(filename, varnames)
%% We already know the file exists by this point.
% Lets check to make sure that file contains all needed variables.

nvars  = length(varnames);
gotone = ones(1,nvars);

for i = 1:nvars
    gotone(i) = nc_var_exists(filename,varnames{i});
    if ( ~ gotone(i) )
        fprintf('\n%s is not a variable in %s\n',varnames{i},filename)
    end
end

if ~ all(gotone)
    error('missing required variable ... exiting')
end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
