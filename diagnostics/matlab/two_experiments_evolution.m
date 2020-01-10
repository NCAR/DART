function two_experiments_evolution(files, titles, obsnames, copy, prpo, varargin)
%% two_experiments_evolution  compares multiple netcdf files created by obs_diag.
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
% USAGE: two_experiments_evolution(files, titles, obsnames, copy, prpo [,varargin])
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
% DateForm : Free-form character string controlling representation of the time axis.
%            See 'help datetick' for discussion and valid values.
%            Example ones are 'mm/dd' and 'dd HH:MM'.
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
% obsnames = {'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_TEMPERATURE'};
%
% copy  = 'spread';   % rmse, spread, totalspread, bias, etc.
% prpo  = 'analysis'; % [analy, analysis, posterior ] == posterior
% prpo  = 'forecast'; % [guess, forecast, prior     ] == prior
%
% two_experiments_evolution(files, titles, obsnames, copy, prpo, 'level', 1)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%%--------------------------------------------------------------------
% Decode,Parse,Check the input
%---------------------------------------------------------------------

default_verbosity  = true;
default_markersize = 12;
default_pause      = false;
default_range      = [NaN NaN];
default_level      = 1;
default_dateform   = 'default';
p = inputParser;

addRequired(p,'files',@iscell);
addRequired(p,'titles',@iscell);
addRequired(p,'obsnames',@iscell);
addRequired(p,'copy',@ischar);
addRequired(p,'prpo',@ischar);

if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'verbose',    default_verbosity,  @islogical);
    addParameter(p,'MarkerSize', default_markersize, @isnumeric);
    addParameter(p,'pause',      default_pause,      @islogical);
    addParameter(p,'range',      default_range,      @isnumeric);
    addParameter(p,'level',      default_level,      @isnumeric);
    addParameter(p,'DateForm',   default_dateform,   @ischar);
else
    addParamValue(p,'verbose',   default_verbosity,  @islogical); %#ok<NVREPL>
    addParamValue(p,'MarkerSize',default_markersize, @isnumeric); %#ok<NVREPL>
    addParamValue(p,'pause',     default_pause,      @islogical); %#ok<NVREPL>
    addParamValue(p,'range',     default_range,      @isnumeric); %#ok<NVREPL>
    addParamValue(p,'level',     default_level,      @isnumeric); %#ok<NVREPL>
    addParamValue(p,'DateForm',  default_dateform,   @ischar);    %#ok<NVREPL>
end

p.parse(files, titles, obsnames, copy, prpo, varargin{:});

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

if (numel(p.Results.range) ~= 2)
    error('range must be an array of length two ... [bottom top]')
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

global figuredata verbose

commondata = check_compatibility(files, prpo, obsnames, copy);
figuredata = set_obsdiag_figure('landscape','numexp',NumExp);
figuredata.MarkerSize = p.Results.MarkerSize;
figuredata.DateForm   = p.Results.DateForm;
verbose    = p.Results.verbose;
plotobj    = cell(NumExp,1);

%%--------------------------------------------------------------------
% Set some static data
%---------------------------------------------------------------------

nvars = length(obsnames);

for ivar = 1:nvars
    fprintf('Working on %s ...\n',obsnames{ivar})
    
    %------------------------------------------------------------------------
    % Plot each region in a separate figure window.
    %------------------------------------------------------------------------
    
    for iregion = 1:commondata.nregions
        
        figure(iregion);
        clf(iregion);
        orient(figuredata.orientation);
        
        %---------------------------------------------------------------------
        % 1) Get the data for each experiment
        % 2) plot the data
        % 3) annotate
        %---------------------------------------------------------------------
        
        for iexp = 1:NumExp
            
            plotobj{iexp} = getvals(files{iexp}, commondata.targets{ivar}, copy, iregion, p.Results.level);
            plotobj{iexp}.title        = titles{iexp};
            plotobj{iexp}.nregions     = commondata.nregions;
            plotobj{iexp}.region_names = commondata.region_names;
            plotobj{iexp}.phase        = commondata.phase;
            
        end
        
        myplot(plotobj);
        
        BottomAnnotation(plotobj)
        
        psfname = sprintf('%s_%s_region%d_ilev%d_evolution_%dexp', ...
            obsnames{ivar}, plotobj{1}.copystring, iregion, p.Results.level, NumExp);
        
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
            targets{i} = sprintf('%s_guess',varnames{i});
            commondata{i}.phase = 'prior';
        case {'analy','analysis','posterior'}
            targets{i} = sprintf('%s_analy',varnames{i});
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
    commondata{i}.ncopies      = nc_dim_info(filenames{i}, 'copy');
    commondata{i}.nobstypes    = nc_dim_info(filenames{i}, 'obstypes');
    commondata{i}.nregions     = nc_dim_info(filenames{i}, 'region');
    commondata{i}.lonlim1      = nc_read_att(filenames{i}, '/','lonlim1');
    commondata{i}.lonlim2      = nc_read_att(filenames{i}, '/','lonlim2');
    commondata{i}.latlim1      = nc_read_att(filenames{i}, '/','latlim1');
    commondata{i}.latlim2      = nc_read_att(filenames{i}, '/','latlim2');
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
    
end

if mystat > 0
    error('The experiments are not compatible ... stopping.')
end

common = commondata{1};


%=====================================================================


function plotdat = getvals(fname, varname, copystring, regionindex, levelindex )
%% Get the data for each experiment
if (exist(fname,'file') ~= 2)
    error('%s does not exist',fname)
end

plotdat         = read_obsdiag_staticdata(fname,copystring);
plotdat.region  = regionindex;
plotdat.varname = varname;

% get appropriate vertical coordinate variable

[dimnames, ~] = nc_var_dims(fname, plotdat.varname);

if ( plotdat.dimensionality == 1 ) % observations on a unit circle, no level
    plotdat.levelindex  = 1;
    plotdat.level       = 1;
    plotdat.level_units = [];
elseif ( strfind(dimnames{2},'surface') > 0 )
    plotdat.levelindex  = 1;
    plotdat.level       = 1;
    plotdat.level_units = 'surface';
    plotdat.level_edges = [];
elseif ( strfind(dimnames{2},'undef') > 0 )
    plotdat.levelindex  = 1;
    plotdat.level       = 1;
    plotdat.level_units = 'undefined';
    plotdat.level_edges = [];
else
    plotdat.levelindex  = levelindex;
    plotdat.level       = ncread(fname, dimnames{2});
    plotdat.level_units = nc_read_att(fname, dimnames{2}, 'units');
    plotdat.level_edges = ncread(fname,sprintf('%s_edges',dimnames{2}));
end

myinfo.diagn_file  = fname;
myinfo.copyindex   = plotdat.copyindex;
myinfo.regionindex = plotdat.region;
myinfo.levelindex  = plotdat.levelindex;
[start, count]     = GetNCindices(myinfo,'diagn',plotdat.varname);
hyperslab          = ncread(fname, plotdat.varname, start, count);
plotdat.data       = squeeze(hyperslab);
plotdat.trusted    = nc_read_att(fname, plotdat.varname, 'TRUSTED');
if (isempty(plotdat.trusted)), plotdat.trusted = 'NO'; end

%% Determine data limits
%  always make sure we have a zero bias line ...

switch copystring
    case {'bias'}
        dmin = min( [ min(plotdat.data) 0.0 ] );
        dmax = max( [ max(plotdat.data) 0.0 ] );
        plotdat.Drange = [ dmin dmax ];
        plotdat.ylabel = sprintf('%s (%s)',copystring, plotdat.biasconv);
    otherwise
        plotdat.Drange = [min(plotdat.data) max(plotdat.data)];
        plotdat.ylabel = copystring;
end

qcvalues = get_qc_values(fname, plotdat.varname, ...
    'regionindex', plotdat.region, ...
    'levelindex',plotdat.levelindex, ...
    'fatal', false, ...
    'verbose', false);

plotdat.nposs         = squeeze(qcvalues.nposs);
plotdat.nused         = squeeze(qcvalues.nused);
plotdat.num_evaluated = squeeze(qcvalues.num_evaluated);

if sum(plotdat.num_evaluated > 0)
    plotdat.assim_eval_string = 'evaluated';
else
    plotdat.assim_eval_string = 'assimilated';
end

%% Set the last of the ranges

plotdat.Nrange = [min(plotdat.nused(:))    max(plotdat.nposs(:))];


%=====================================================================


function myplot(plotobj)
%% myplot Creates a graphic for one region

global figuredata

Nexp    = length(plotobj);

%% Create the background

ax1   = subplot('position',figuredata.position);
set(ax1,'YAxisLocation','left','FontSize',figuredata.fontsize)

%% draw the results of the experiments, priors and posteriors
%  each with their own line type.
hd     = [];   % handle to an unknown number of data lines
legstr = {[]}; % strings for the legend

for i = 1:Nexp
    hd(i) = line(plotobj{i}.bincenters, plotobj{i}.data, 'Parent', ax1); %#ok<AGROW>
    
    set(hd(i), 'Color',    figuredata.expcolors{i}, ...
        'Marker',          figuredata.expsymbols{i}, ...
        'MarkerFaceColor', figuredata.expcolors{i}, ...
        'MarkerSize',      figuredata.MarkerSize, ...
        'LineStyle',       figuredata.prpolines{1}, ...
        'LineWidth',       figuredata.linewidth);
    
    % calculate the weighted mean for a summary. Each experiment
    % may use different numbers of observations, so a simple mean
    % may be misleading.
    
    N = sum(plotobj{i}.nused,'omitnan');
    X = sum(plotobj{i}.data .* plotobj{i}.nused,'omitnan')/N;
    
    legstr{i} = sprintf('%s ... mean = %s',plotobj{i}.title,num2str(X));
end

% Plot a bias line.

switch plotobj{1}.copystring
    case {'bias'}
        zeroline = line(get(ax1,'XLim'),[0 0],'Color',[200 200 200]/255,'Parent',ax1);
        set(zeroline,'LineWidth',2.5,'LineStyle','-')
    otherwise
end

% effort to use user-supplied value for time labelling or
% make a stab at a useful default.

set_time_axis('x', plotobj{i}.bincenters, figuredata.DateForm);

% Create another axes to use for plotting the observation counts
% using a black axis because there is no single observation color.

ax2 = axes( ...
    'Position',get(ax1,'Position'), ...
    'FontSize',get(ax1,'FontSize'), ...
    'XColor'  ,get(ax1,'XColor'), ...
    'XLim'    ,get(ax1,'XLim'), ...
    'XTick'   ,get(ax1,'XTick'), ...
    'YDir'    ,get(ax1,'YDir'), ...
    'Color'   ,'none', ...
    'YColor'  ,'k', ...
    'XAxisLocation','top', ...
    'YAxisLocation','right');

% Plot the data, which sets the range of the axis
for i = 1:Nexp
    ax2h1 = line(plotobj{i}.bincenters, plotobj{i}.nposs, 'Parent',ax2);
    ax2h2 = line(plotobj{i}.bincenters, plotobj{i}.nused, 'Parent',ax2);
    
    set(ax2h1,'LineStyle','none', ...
        'Color',      figuredata.expcolors{i}, ...
        'Marker',     figuredata.obs_marker, ...
        'MarkerSize', figuredata.MarkerSize);
    
    set(ax2h2,'LineStyle','none', ...
        'Color',      figuredata.expcolors{i}, ...
        'Marker',     figuredata.ges_marker, ...
        'MarkerSize', figuredata.MarkerSize);
end

% turn off topside X tick labels (clashes with title)
% use the same Y ticks, but find the right label values
set(ax2, 'XTicklabel', []);
matchingYticks(ax1,ax2);

% Annotate. Trying to maximize content, minimize clutter.
annotate( ax1, ax2, plotobj{1})

lh = legend(hd,legstr);
set(lh,'Interpreter','none','Box','off','FontSize',figuredata.fontsize);

if verLessThan('matlab','R2017a')
    % Convince Matlab to not autoupdate the legend with each new line.
    % Before 2017a, this was the default behavior, so do nothing.
    % We do not want to add the bias line to the legend, for example.
else
    lh.AutoUpdate = 'off';
end

%=====================================================================


function annotate(ax1, ax2, plotobj)
%% One figure ... everything gets annotated.

global figuredata

set(get(ax1,'Xlabel'),'String',plotobj.timespan, ...
    'Interpreter','none','FontSize',figuredata.fontsize)

ylabel = sprintf('%s (%s)',plotobj.ylabel, plotobj.phase);

set(get(ax1,'Ylabel'),'String',ylabel, ...
    'Interpreter','none','FontSize',figuredata.fontsize)

string1 = sprintf('# of obs (o=possible, %s=%s)', '\ast', plotobj.assim_eval_string);

set(get(ax2,'Ylabel'), 'String', string1, 'FontSize', figuredata.fontsize)

if ( isempty(plotobj.level_units) )
    th = title({deblank(plotobj.region_names(plotobj.region,:)), ...
        sprintf('%s @ %f', plotobj.varname, plotobj.level(plotobj.levelindex))});
else
    th = title({deblank(plotobj.region_names(plotobj.region,:)), ...
        sprintf('%s @ %f %s', plotobj.varname, plotobj.level(plotobj.levelindex), ...
        plotobj.level_units)});
end
set(th,'Interpreter','none','FontSize',figuredata.fontsize,'FontWeight','bold');


%=====================================================================


function BottomAnnotation(plotstruct)
%% annotates the filename containing the data being plotted

Nexp    = length(plotstruct);
dy      = 1.0/(Nexp+2);
yheight = 0.0225*(Nexp+1);

subplot('position',[0.10 0.01 0.8 yheight])
axis off

for ifile = 1:Nexp
    main = plotstruct{ifile}.fname;
    
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
    
    ty = 1.0 - (ifile+1)*dy;
    h = text(0.0, ty, string1);
    set(h, 'Interpreter', 'none', ...
        'HorizontalAlignment','left', ...
        'FontSize', 8);
    
    % If the observation is trusted for this experiment, annotate as such.
    
    switch lower(plotstruct{ifile}.trusted)
        case 'true'
            ty = 1.0 - Nexp*dy + ifile*dy*2;
            h = text(0.0,ty,sprintf('TRUSTED OBSERVATION in %s',plotstruct{ifile}.title));
            set(h,'FontSize',20,'Interpreter','none','HorizontalAlignment','left')
        otherwise
    end
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
