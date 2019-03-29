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
% USAGE: plotdat = plot_profile(fname, copy);
%
% fname    :  netcdf file produced by 'obs_diag'
%
% copy     : string defining the metric of interest. 'rmse', 'spread', etc.
%            Possible values are available in the netcdf 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% obsname  : Optional. If present, The strings of each observation type to plot.
%            Each observation type will be plotted in a separate graphic.
%            Default is to plot all available observation types.
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

default_obsname = 'none';
default_verbosity = 'yes';
p = inputParser;

addRequired(p,'fname',@ischar);
addRequired(p,'copy',@ischar);
if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'obsname',default_obsname,@ischar);
    addParameter(p,'verbose',default_verbosity,@ischar);
else
    addParamValue(p,'obsname',default_obsname,@ischar);
    addParamValue(p,'verbose',default_verbosity,@ischar);
end

p.parse(fname, copy, varargin{:});

% if you want to echo the input
% disp(['fname   : ', p.Results.fname])
% disp(['copy    : ', p.Results.copy])
% disp(['obsname : ', p.Results.obsname])

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

global verbose
if (strncmpi(p.Results.verbose,'y',1))
    verbose = 1;
else
    verbose = 0;
end

if (exist(fname,'file') ~= 2)
    error('file/fname <%s> does not exist',fname)
end

%%--------------------------------------------------------------------
% Harvest plotting info/metadata from netcdf file.
%---------------------------------------------------------------------

plotdat.fname         = fname;
plotdat.copystring    = copy;

plotdat.binseparation = nc_read_att(fname, '/', 'bin_separation');
plotdat.binwidth      = nc_read_att(fname, '/', 'bin_width');
time_to_skip          = nc_read_att(fname, '/', 'time_to_skip');
plotdat.lonlim1       = nc_read_att(fname, '/', 'lonlim1');
plotdat.lonlim2       = nc_read_att(fname, '/', 'lonlim2');
plotdat.latlim1       = nc_read_att(fname, '/', 'latlim1');
plotdat.latlim2       = nc_read_att(fname, '/', 'latlim2');
plotdat.biasconv      = nc_read_att(fname, '/', 'bias_convention');

plotdat.mlevel        = local_ncread(fname, 'mlevel');
plotdat.plevel        = local_ncread(fname, 'plevel');
plotdat.plevel_edges  = local_ncread(fname, 'plevel_edges');
plotdat.hlevel        = local_ncread(fname, 'hlevel');
plotdat.hlevel_edges  = local_ncread(fname, 'hlevel_edges');
plotdat.bincenters    = ncread(fname, 'time');
plotdat.binedges      = ncread(fname, 'time_bounds');
plotdat.region_names  = strtrim(ncread(fname, 'region_names')');
[plotdat.nregions,~]  = nc_dim_info(fname,'region');

% Coordinate between time types and dates

timeunits             = nc_read_att(fname,'time','units');
timebase              = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin            = datenum(timebase(1),timebase(2),timebase(3));
timefloats            = zeros(size(time_to_skip));  % stupid int32 type conversion
timefloats(:)         = time_to_skip(:);
skip_seconds          = timefloats(4)*3600 + timefloats(5)*60 + timefloats(6);
iskip                 = timefloats(3) + skip_seconds/86400.0;

% Set up a structure to use for plotting

plotdat.bincenters    = plotdat.bincenters + timeorigin;
plotdat.binedges      = plotdat.binedges   + timeorigin;
plotdat.Nbins         = length(plotdat.bincenters);
plotdat.toff          = plotdat.binedges(1) + iskip;
plotdat.timespan      = sprintf('%s through %s', datestr(plotdat.toff), ...
    datestr(max(plotdat.binedges(:))));
plotdat.xlabel        = sprintf('%s',copy);

[plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(fname);
[plotdat.varnames,    plotdat.vardims]    = FindVerticalVars(plotdat);

plotdat.nvars         = length(plotdat.varnames);
plotdat.copyindex     = get_copy_index(fname,copy);
plotdat.rmseindex     = get_copy_index(fname,'rmse');
plotdat.Npossindex    = get_copy_index(fname,'Nposs');
plotdat.Nusedindex    = get_copy_index(fname,'Nused');
plotdat.NQC0index     = get_copy_index(fname,'N_DARTqc_0');
plotdat.NQC1index     = get_copy_index(fname,'N_DARTqc_1');
plotdat.NQC2index     = get_copy_index(fname,'N_DARTqc_2');
plotdat.NQC3index     = get_copy_index(fname,'N_DARTqc_3');
plotdat.NQC4index     = get_copy_index(fname,'N_DARTqc_4');
plotdat.NQC5index     = get_copy_index(fname,'N_DARTqc_5');
plotdat.NQC6index     = get_copy_index(fname,'N_DARTqc_6');
plotdat.NQC7index     = get_copy_index(fname,'N_DARTqc_7');
plotdat.NQC8index     = get_copy_index(fname,'N_DARTqc_8','fatal',false);

global figuredata
figuredata = setfigure();

%%---------------------------------------------------------------------
% Loop around (copy-level-region) observation types
%----------------------------------------------------------------------

% Either use all the variables or just the one optionally specified.

if strcmp(p.Results.obsname,'none')
    varlist = 1:plotdat.nvars;
else
    varlist = find (strcmpi(p.Results.obsname,plotdat.varnames));
    if isempty(varlist)
        error('%s is not in the list of observations',p.Results.obsname)
    end
end

for ivar = varlist

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

    [level_org, level_units, nlevels, level_edges, Yrange] = FindVerticalInfo(fname, plotdat.guessvar);
    plotdat.level_org   = level_org;
    plotdat.level_units = level_units;
    plotdat.nlevels     = nlevels;
    plotdat.level_edges = level_edges;
    plotdat.Yrange      = Yrange;

    % Matlab likes strictly ASCENDING order for things to be plotted,
    % then you can impose the direction. The data is stored in the original
    % order, so the sort indices are saved to reorder the data.

    if (plotdat.level_org(1) > plotdat.level_org(plotdat.nlevels))
        plotdat.YDir = 'reverse';
    else
        plotdat.YDir = 'normal';
    end

    % Add error-checking for output from older versions of obs_diag.

    [levels, indices]   = sort(plotdat.level_org);
    plotdat.level       = unique(levels);
    if (length(plotdat.level) ~= length(levels))
        error('There is a duplicated value in the array specifying the levels - must change your input.nml and rerun obs_diag')
    end

    plotdat.indices     = indices;
    level_edges         = sort(plotdat.level_edges);
    plotdat.level_edges = level_edges;

    % The rest of this script was written for the third-party netcdf
    % support. Matlab's native ncread transposes the variables, so I have to
    % permute them back to the expected storage order.

    guess = ncread(fname, plotdat.guessvar);
    analy = local_ncread(fname, plotdat.analyvar);
    if (isempty(analy))
        plotdat.has_analysis = 0;
        analy = guess;  % make the variable the same shape as guess
        analy(:) = NaN; % and fill it with nothing
    else
        plotdat.has_analysis = 1;
    end
    rank  = length(size(guess));
    guess = permute(guess,rank:-1:1);
    analy = permute(analy,rank:-1:1);

    % singleton dimensions are auto-squeezed - which is unfortunate.
    % We want these things to be 3D. [copy-level-region]

    if ( plotdat.nlevels == 1 )
        bob(:,1,:) = guess;
        ted(:,1,:) = analy;
        guess = bob; clear bob
        analy = ted; clear ted
    end

    % check to see if there is anything to plot
    % The number possible is decreased by the number of observations
    % rejected by namelist control.

    % todo write a function to print this information given a structure
    % with the indices. Many of these are not needed, but a report is good

    if (verbose)
        fprintf('\n')
        fprintf('%10d %s observations had DART QC of 4 (all regions).\n', ...
            sum(sum(guess(plotdat.NQC4index, :,:))),plotdat.myvarname)
        fprintf('%10d %s observations had DART QC of 5 (all regions).\n', ...
            sum(sum(guess(plotdat.NQC5index, :,:))),plotdat.myvarname)
        fprintf('%10d %s observations had DART QC of 6 (all regions).\n', ...
            sum(sum(guess(plotdat.NQC6index, :,:))),plotdat.myvarname)
        fprintf('%10d %s observations had DART QC of 7 (all regions).\n', ...
            sum(sum(guess(plotdat.NQC7index, :,:))),plotdat.myvarname)
        if (plotdat.NQC8index > 0)
            fprintf('%10d %s observations had DART QC of 8 (all regions).\n', ...
                sum(sum(guess(plotdat.NQC8index, :,:))),plotdat.myvarname)
        end
        fprintf('\n')
    end

    nposs = sum(guess(plotdat.Npossindex,:,:)) - ...
        sum(guess(plotdat.NQC5index ,:,:)) - ...
        sum(guess(plotdat.NQC6index ,:,:));

    if ( sum(nposs(:)) < 1 )
        fprintf('No obs for %s...  skipping\n', plotdat.varnames{ivar})
        continue
    end

    plotdat.ges_copy   = guess(plotdat.copyindex,  :, :);
    plotdat.ges_rmse   = guess(plotdat.rmseindex,  :, :);
    plotdat.ges_Nqc0   = guess(plotdat.NQC0index,  :, :);
    plotdat.ges_Nqc1   = guess(plotdat.NQC1index,  :, :);
    plotdat.ges_Nqc2   = guess(plotdat.NQC2index,  :, :);
    plotdat.ges_Nqc3   = guess(plotdat.NQC3index,  :, :);
    plotdat.ges_Nqc4   = guess(plotdat.NQC4index,  :, :);
    plotdat.ges_Nqc5   = guess(plotdat.NQC5index,  :, :);
    plotdat.ges_Nqc6   = guess(plotdat.NQC6index,  :, :);
    plotdat.ges_Nqc7   = guess(plotdat.NQC7index,  :, :);
    if (plotdat.NQC8index > 0)
        plotdat.ges_Nqc8   = guess(plotdat.NQC8index,  :, :);
    end

    plotdat.ges_Nused  = guess(plotdat.Nusedindex, :, :);
    plotdat.ges_Nposs  = guess(plotdat.Npossindex, :, :) - ...
        plotdat.ges_Nqc5 - plotdat.ges_Nqc6;

    plotdat.anl_copy   = analy(plotdat.copyindex,  :, :);
    plotdat.anl_rmse   = analy(plotdat.rmseindex,  :, :);
    plotdat.anl_Nqc0   = analy(plotdat.NQC0index,  :, :);
    plotdat.anl_Nqc1   = analy(plotdat.NQC1index,  :, :);
    plotdat.anl_Nqc2   = analy(plotdat.NQC2index,  :, :);
    plotdat.anl_Nqc3   = analy(plotdat.NQC3index,  :, :);
    plotdat.anl_Nqc4   = analy(plotdat.NQC4index,  :, :);
    plotdat.anl_Nqc5   = analy(plotdat.NQC5index,  :, :);
    plotdat.anl_Nqc6   = analy(plotdat.NQC6index,  :, :);
    plotdat.anl_Nqc7   = analy(plotdat.NQC7index,  :, :);
    if (plotdat.NQC8index > 0)
        plotdat.anl_Nqc8   = analy(plotdat.NQC8index,  :, :);
    end

    plotdat.anl_Nused  = analy(plotdat.Nusedindex, :, :);
    plotdat.anl_Nposs  = analy(plotdat.Npossindex, :, :) - ...
        plotdat.anl_Nqc5 - plotdat.anl_Nqc6;

    plotdat.Xrange     = FindRange(plotdat);

    % plot by region - each in its own figure.

    for iregion = 1:plotdat.nregions
        figure(iregion); clf(iregion); orient(figuredata.orientation);
        plotdat.region   = iregion;
        plotdat.myregion = deblank(plotdat.region_names(iregion,:));
        myplot(plotdat);
        BottomAnnotation(fname)

        psfname = sprintf('%s_%s_profile_region%d', plotdat.varnames{ivar}, ...
            plotdat.copystring, iregion);
        print(gcf,'-dpdf',psfname);
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

% add type of vertical coordinate info for adjusting axes to accomodate legend

Stripes(plotdat.Xrange, plotdat.level_edges, plotdat.level_units);
set(ax1, 'YDir', plotdat.YDir, 'YTick', plotdat.level, 'Layer', 'top')
set(ax1,'YAxisLocation','left','FontSize',figuredata.fontsize)

% draw the result of the experiment

[hprior, prior_legstr] = plot_quantity('prior', plotdat);

nobs_poss       =     plotdat.ges_Nposs(:,:,plotdat.region);
nobs_used_prior =     plotdat.ges_Nused(:,:,plotdat.region);
nobs_used_poste =     plotdat.anl_Nused(:,:,plotdat.region);
anl_Nposs       = sum(plotdat.anl_Nposs(:,:,plotdat.region));

if (verbose)
    fprintf('region %d nobs_poss %d prior %d poste %d\n', plotdat.region, ...
        sum(nobs_poss), sum(nobs_used_prior), sum(nobs_used_poste))
    fprintf('region %d prior_legstr %s\n\n', plotdat.region, prior_legstr)
end

if ( isfinite(anl_Nposs) )
    [hposte, poste_legstr] = plot_quantity('posterior', plotdat);
    h = legend([hprior, hposte], prior_legstr, poste_legstr);
else
    h = legend(hprior,prior_legstr);
end

switch plotdat.copystring
    case {'bias'}
        zeroline = line([0 0],plotdat.Yrange,'Color',[0 100 0]/255,'Parent',ax1);
        set(zeroline,'LineWidth',2.5,'LineStyle','-')
        plotdat.xlabel = sprintf('%s (%s)',plotdat.copystring, plotdat.biasconv);
    otherwise
end

set(h,'Interpreter','none','Box','off','Location','NorthWest')

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

ax2h1 = line(nobs_poss      ,plotdat.level,'Color',figuredata.obs_color,'Parent',ax2);
ax2h2 = line(nobs_used_prior,plotdat.level,'Color',figuredata.ges_color,'Parent',ax2);
set(ax2h1,'LineStyle','none','Marker','o');
set(ax2h2,'LineStyle','none','Marker',figuredata.ges_marker);

if (isfinite(sum(anl_Nposs)))
    ax2h3 = line(nobs_used_poste,plotdat.level,'Color',figuredata.anl_color,'Parent',ax2);
    set(ax2h3,'LineStyle','none','Marker',figuredata.anl_marker);
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

nevaluated = sum(plotdat.ges_Nqc1(:) + plotdat.ges_Nqc3(:));

if (nevaluated > 0)
    string1 = sprintf('# of obs (o=possible; %s=evaluated) x %d', ...
        '\ast,\diamondsuit', uint32(xscale));
else
    string1 = sprintf('# of obs (o=possible; %s=assimilated) x %d', ...
        '\ast,\diamondsuit', uint32(xscale));
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

if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
    error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

global verbose

j = 0;
basenames = struct([]);
basedims  = struct([]);

for i = 1:length(x.allvarnames)
    dimnames = lower(x.allvardims{i});
    if (isempty(strfind(dimnames,'time')))
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


function [level_org, level_units, nlevels, level_edges, Yrange] = FindVerticalInfo(fname,varname)
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

level_org   = ncread(   fname,varinfo.Dimensions(leveldim).Name);
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


function figdata = setfigure()
%%
%  figure out a page layout
%  extra space at the bottom for the date/file annotation
%  extra space at the top because the titles have multiple lines

orientation   = 'tall';
position      = [0.15 0.12 0.7 0.75];
fontsize      = 16;
linewidth     = 2.5;
markersize    = 8.0;
ges_color     = [  0/255 128/255   0/255]; % prior_green
anl_color     = [  0/255   0/255 255/255]; % poste_blue
obs_color     = [215/255  10/255  83/255]; % obs_red
ges_marker    = '*';
anl_marker    = 'd';
ges_linestyle = '-';
anl_linestyle = '-';

figdata = struct( ...
    'expcolors',  {{'k','r','b','m','g','c','y'}}, ...
    'expsymbols', {{'o','s','d','p','h','s','*'}}, ...
    'prpolines',  {{'-','--'}},  ...
    'position'     , position, ...
    'fontsize'     , fontsize,  ...
    'orientation'  , orientation, ...
    'linewidth'    , linewidth, ...
    'markersize'   , markersize, ...
    'ges_color'    , ges_color, ...
    'anl_color'    , anl_color, ...
    'obs_color'    , obs_color, ...
    'ges_marker'   , ges_marker, ...
    'anl_marker'   , anl_marker, ...
    'ges_linestyle', ges_linestyle, ...
    'anl_linestyle', anl_linestyle );


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
        data      = plotdat.ges_copy( :,:,plotdat.region);
        Nused     = plotdat.ges_Nused(:,:,plotdat.region);
        color     = figuredata.ges_color;
        marker    = figuredata.ges_marker;
        linestyle = figuredata.ges_linestyle;
        linewidth = figuredata.linewidth;
        string1   = 'forecast:';
    case 'posterior'
        data      = plotdat.anl_copy( :,:,plotdat.region);
        Nused     = plotdat.anl_Nused(:,:,plotdat.region);
        color     = figuredata.anl_color;
        marker    = figuredata.anl_marker;
        linestyle = figuredata.anl_linestyle;
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

h = line(data, plotdat.level);
set(h, 'LineStyle',  linestyle, ...
    'LineWidth',  linewidth, ...
    'Color',      color, ...
    'Marker',     marker,    ...
    'MarkerSize', 3*linewidth);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
