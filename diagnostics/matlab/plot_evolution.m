function plotdat = plot_evolution(fname, copy, varargin)
%% plot_evolution plots the temporal evolution of the observation-space quantities for all possible levels, all possible variables.
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
% USAGE: plotdat = plot_evolution(fname, copy);
%
% fname    :  netcdf file produced by 'obs_diag'
%
% copy     : string defining the metric of interest. 'rmse', 'spread', etc.
%            Possible values are available in the netcdf 'CopyMetaData' variable.
%            (ncdump -v CopyMetaData obs_diag_output.nc)
%
% varargin: optional, parameter-value pairs. Supported parameters are described below.
%
% obsname  : The strings of each observation type to plot.
%            Each observation type will be plotted in a separate graphic.
%            Default is to plot all available observation types.
%
% level    :  'level' index. Default is to plot all levels.
%
% range    : 'range' of the value being plotted. Default is to
%                automatically determine range based on the data values.
%
% OUTPUT: 'plotdat' is a structure containing what was last plotted.
%         A postscript file containing a page for each level - each region.
%         The other file is a simple text file containing summary information
%         about how many observations were assimilated, how many were available, etc.
%         Both of these filenames contain the observation type,
%         copy and region as part of the name.
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
% plotdat    = plot_evolution(fname, copy, 'obsname', 'RADIOSONDE_TEMPERATURE', ...
%                             'level', 4, 'range', [0 10]);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

default_obsname = 'none';
default_range = [NaN NaN];
default_level = -1;
default_verbosity = 'yes';
default_markersize = 8;
default_pause = 'no';
p = inputParser;

addRequired(p,'fname',@ischar);
addRequired(p,'copy',@ischar);
if (exist('inputParser/addParameter','file') == 2)
    addParameter(p,'obsname',    default_obsname,    @ischar);
    addParameter(p,'range',      default_range,      @isnumeric);
    addParameter(p,'level',      default_level,      @isnumeric);
    addParameter(p,'verbose',    default_verbosity,  @ischar);
    addParameter(p,'MarkerSize', default_markersize, @isnumeric);
    addParameter(p,'pause',      default_pause,      @ischar);
else
    addParamValue(p,'obsname',   default_obsname,    @ischar);
    addParamValue(p,'range',     default_range,      @isnumeric);
    addParamValue(p,'level',     default_level,      @isnumeric);
    addParamValue(p,'verbose',   default_verbosity,  @ischar);
    addParamValue(p,'MarkerSize',default_markersize, @isnumeric);
    addParamValue(p,'pause',     default_pause,      @ischar);
end
p.parse(fname, copy, varargin{:});

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

plotdat.fname         = fname;
plotdat.copystring    = copy;
plotdat.bincenters    = ncread(fname,'time');
plotdat.binedges      = ncread(fname,'time_bounds');
plotdat.mlevel        = local_ncread(fname,'mlevel');
plotdat.plevel        = local_ncread(fname,'plevel');
plotdat.plevel_edges  = local_ncread(fname,'plevel_edges');
plotdat.hlevel        = local_ncread(fname,'hlevel');
plotdat.hlevel_edges  = local_ncread(fname,'hlevel_edges');
[plotdat.ncopies, ~]  = nc_dim_info(fname,'copy');
[plotdat.nregions, ~] = nc_dim_info(fname,'region');
plotdat.region_names  = strtrim(ncread(fname,'region_names')');

dimensionality        = nc_read_att(fname, '/', 'LocationRank');
plotdat.binseparation = nc_read_att(fname, '/', 'bin_separation');
plotdat.binwidth      = nc_read_att(fname, '/', 'bin_width');
time_to_skip          = nc_read_att(fname, '/', 'time_to_skip');
plotdat.lonlim1       = nc_read_att(fname, '/', 'lonlim1');
plotdat.lonlim2       = nc_read_att(fname, '/', 'lonlim2');
plotdat.latlim1       = nc_read_att(fname, '/', 'latlim1');
plotdat.latlim2       = nc_read_att(fname, '/', 'latlim2');
plotdat.biasconv      = nc_read_att(fname, '/', 'bias_convention');

% Coordinate between time types and dates

%calendar     = nc_read_att(fname,'time','calendar');
timeunits    = nc_read_att(fname,'time','units');
timebase     = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin   = datenum(timebase(1),timebase(2),timebase(3));
if ( isempty(time_to_skip) == 1)
    iskip = 0;
elseif ( numel(time_to_skip) == 6)
    skip_seconds = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6);
    iskip        = time_to_skip(3) + skip_seconds/86400;
else
    error('time_to_skip variable has unusual length. Should be either 0 or 6.')
end

% set up a structure with all static plotting components

plotdat.bincenters = plotdat.bincenters + timeorigin;
plotdat.binedges   = plotdat.binedges   + timeorigin;
plotdat.Nbins      = length(plotdat.bincenters);
plotdat.toff       = plotdat.bincenters(1) + iskip;

if (nvars == 0)
    [plotdat.allvarnames, plotdat.allvardims] = get_varsNdims(fname);
    [plotdat.varnames,    plotdat.vardims]    = FindTemporalVars(plotdat);
    plotdat.nvars       = length(plotdat.varnames);
else
    plotdat.varnames{1} = obsname;
    plotdat.nvars       = nvars;
end

plotdat.copyindex   = get_copy_index(fname,copy);
plotdat.Npossindex  = get_copy_index(fname,'Nposs');
plotdat.Nusedindex  = get_copy_index(fname,'Nused');
plotdat.NQC0index   = get_copy_index(fname,'N_DARTqc_0');
plotdat.NQC1index   = get_copy_index(fname,'N_DARTqc_1');
plotdat.NQC2index   = get_copy_index(fname,'N_DARTqc_2');
plotdat.NQC3index   = get_copy_index(fname,'N_DARTqc_3');
plotdat.NQC4index   = get_copy_index(fname,'N_DARTqc_4');
plotdat.NQC5index   = get_copy_index(fname,'N_DARTqc_5');
plotdat.NQC6index   = get_copy_index(fname,'N_DARTqc_6');
plotdat.NQC7index   = get_copy_index(fname,'N_DARTqc_7');
plotdat.NQC8index   = get_copy_index(fname,'N_DARTqc_8','fatal',false);

global figuredata
figuredata = setfigure();
figuredata.MarkerSize = p.Results.MarkerSize;

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
        fprintf('Removing %s from the current directory.\n',psfname{iregion})
        system(sprintf('rm %s',psfname{iregion}));
    end

    % remove any existing log file -

    lgfname = sprintf('%s_%s_obscount.txt',plotdat.varnames{ivar},plotdat.copystring);
    fprintf('Removing %s from the current directory.\n',lgfname)
    system(sprintf('rm %s',lgfname));
    logfid = fopen(lgfname,'wt');
    fprintf(logfid,'%s\n',lgfname);

    % get appropriate vertical coordinate variable

    [dimnames, ~] = nc_var_dims(fname, plotdat.guessvar);

    if ( dimensionality == 1 ) % observations on a unit circle, no level
        plotdat.level = 1;
        plotdat.level_units = [];
    elseif ( strfind(dimnames{2},'surface') > 0 )
        plotdat.level       = 1;
        plotdat.level_units = 'surface';
    elseif ( strfind(dimnames{2},'undef') > 0 )
        plotdat.level       = 1;
        plotdat.level_units = 'undefined';
    else
        plotdat.level       = ncread(fname, dimnames{2});
        plotdat.level_units = nc_read_att(fname, dimnames{2}, 'units');
    end
    plotdat.nlevels = length(plotdat.level);

    % Here is the tricky part. Singleton dimensions are auto-squeezed ...
    % single levels, single regions ...

    guess_raw = ncread(fname, plotdat.guessvar);
    guess_raw = permute(guess_raw,length(size(guess_raw)):-1:1);
    guess = reshape(guess_raw, plotdat.Nbins,   plotdat.ncopies, ...
        plotdat.nlevels, plotdat.nregions);

    analy_raw = local_ncread(fname, plotdat.analyvar);
    if ( isempty(analy_raw) )
        % force analysis to be the same shape as the guess,
        % and full of NaNs instead of 0s.
        analy = guess;
        analy(:) = NaN;
    else
        analy_raw = permute(analy_raw,length(size(analy_raw)):-1:1);
        analy = reshape(analy_raw, plotdat.Nbins,   plotdat.ncopies, ...
            plotdat.nlevels, plotdat.nregions);
    end

    % check to see if there is anything to plot
    % The number possible is decreased by the number of observations
    % rejected by namelist control.

    qc4  = guess(:,plotdat.NQC4index ,:,:);
    qc5  = guess(:,plotdat.NQC5index ,:,:);
    qc6  = guess(:,plotdat.NQC6index ,:,:);
    qc7  = guess(:,plotdat.NQC7index ,:,:);
    qnp  = guess(:,plotdat.Npossindex,:,:);

    nqc4  = sum(qc4(:));
    nqc5  = sum(qc5(:));
    nqc6  = sum(qc6(:));
    nqc7  = sum(qc7(:));
    nposs = sum(qnp(:)) - nqc5 - nqc6;

    clear qc4 qc5 qc6 qc7 qnp

    if (verbose)
        fprintf('\n')
        fprintf('%10d %s observations had DART QC of 5 (all levels, all regions).\n', ...
            nqc5, plotdat.myvarname)
        fprintf('%10d %s observations had DART QC of 6 (all levels, all regions).\n', ...
            nqc6, plotdat.myvarname)
        fprintf('%10d %s observations remain  possible (all levels, all regions).\n', ...
            nposs, plotdat.myvarname)
        fprintf('%10d %s observations had DART QC of 4 (all levels, all regions).\n', ...
            nqc4, plotdat.myvarname)
        fprintf('%10d %s observations had DART QC of 7 (all levels, all regions).\n', ...
            nqc7, plotdat.myvarname)
        fprintf('\n')
    end

    if ( sum(nposs(:)) < 1 )
        fprintf('no obs for %s...  skipping\n', plotdat.varnames{ivar})
        continue
    end

    if (p.Results.level < 0)
        wantedlevels = 1:plotdat.nlevels;
    else
        wantedlevels = p.Results.level;
    end

    for ilevel = wantedlevels

        plotdat.mylevel = ilevel;

        % summarize the observation counts in the log file

        fprintf(logfid,'\nlevel %d %f %s\n',ilevel,plotdat.level(ilevel),plotdat.level_units);
        plotdat.ges_Nqc0  = guess(:,plotdat.NQC0index  ,ilevel,:);
        plotdat.anl_Nqc0  = analy(:,plotdat.NQC0index  ,ilevel,:);
        fprintf(logfid,'DART QC == 0, prior/post %d %d\n',sum(plotdat.ges_Nqc0(:)), ...
            sum(plotdat.anl_Nqc0(:)));

        plotdat.ges_Nqc1  = guess(:,plotdat.NQC1index  ,ilevel,:);
        plotdat.anl_Nqc1  = analy(:,plotdat.NQC1index  ,ilevel,:);
        fprintf(logfid,'DART QC == 1, prior/post %d %d\n',sum(plotdat.ges_Nqc1(:)), ...
            sum(plotdat.anl_Nqc1(:)));

        plotdat.ges_Nqc2  = guess(:,plotdat.NQC2index  ,ilevel,:);
        plotdat.anl_Nqc2  = analy(:,plotdat.NQC2index  ,ilevel,:);
        fprintf(logfid,'DART QC == 2, prior/post %d %d\n',sum(plotdat.ges_Nqc2(:)), ...
            sum(plotdat.anl_Nqc2(:)));

        plotdat.ges_Nqc3  = guess(:,plotdat.NQC3index  ,ilevel,:);
        plotdat.anl_Nqc3  = analy(:,plotdat.NQC3index  ,ilevel,:);
        fprintf(logfid,'DART QC == 3, prior/post %d %d\n',sum(plotdat.ges_Nqc3(:)), ...
            sum(plotdat.anl_Nqc3(:)));

        plotdat.ges_Nqc4  = guess(:,plotdat.NQC4index  ,ilevel,:);
        plotdat.anl_Nqc4  = analy(:,plotdat.NQC4index  ,ilevel,:);
        fprintf(logfid,'DART QC == 4, prior/post %d %d\n',sum(plotdat.ges_Nqc4(:)), ...
            sum(plotdat.anl_Nqc4(:)));

        plotdat.ges_Nqc5  = guess(:,plotdat.NQC5index  ,ilevel,:);
        plotdat.anl_Nqc5  = analy(:,plotdat.NQC5index  ,ilevel,:);
        fprintf(logfid,'DART QC == 5, prior/post %d %d\n',sum(plotdat.ges_Nqc5(:)), ...
            sum(plotdat.anl_Nqc5(:)));

        plotdat.ges_Nqc6  = guess(:,plotdat.NQC6index  ,ilevel,:);
        plotdat.anl_Nqc6  = analy(:,plotdat.NQC6index  ,ilevel,:);
        fprintf(logfid,'DART QC == 6, prior/post %d %d\n',sum(plotdat.ges_Nqc6(:)), ...
            sum(plotdat.anl_Nqc6(:)));

        plotdat.ges_Nqc7  = guess(:,plotdat.NQC7index  ,ilevel,:);
        plotdat.anl_Nqc7  = analy(:,plotdat.NQC7index  ,ilevel,:);
        fprintf(logfid,'DART QC == 7, prior/post %d %d\n',sum(plotdat.ges_Nqc7(:)), ...
            sum(plotdat.anl_Nqc7(:)));

        if (plotdat.NQC8index > 0)
            plotdat.ges_Nqc8  = guess(:,plotdat.NQC8index  ,ilevel,:);
            plotdat.anl_Nqc8  = analy(:,plotdat.NQC8index  ,ilevel,:);
            fprintf(logfid,'DART QC == 8, prior/post %d %d\n',sum(plotdat.ges_Nqc8(:)), ...
                sum(plotdat.anl_Nqc8(:)));
        end

        plotdat.ges_Nposs = guess(:,plotdat.Npossindex, ilevel,:) - ...
            plotdat.ges_Nqc5 - plotdat.ges_Nqc6;
        plotdat.anl_Nposs = analy(:,plotdat.Npossindex, ilevel,:) - ...
            plotdat.anl_Nqc5 - plotdat.anl_Nqc6;
        fprintf(logfid,'# obs poss,   prior/post %d %d\n',sum(plotdat.ges_Nposs(:)), ...
            sum(plotdat.anl_Nposs(:)));

        plotdat.ges_Nused = guess(:,plotdat.Nusedindex, ilevel,:);
        plotdat.anl_Nused = analy(:,plotdat.Nusedindex, ilevel,:);
        fprintf(logfid,'# obs used,   prior/post %d %d\n',sum(plotdat.ges_Nused(:)), ...
            sum(plotdat.anl_Nused(:)));

        plotdat.ges_copy  = guess(:,plotdat.copyindex,  ilevel,:);
        plotdat.anl_copy  = analy(:,plotdat.copyindex,  ilevel,:);

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
                plotdat.title    = plotdat.myvarname;
            else
                plotdat.title    = sprintf('%s @ %d %s',    ...
                    plotdat.myvarname,     ...
                    plotdat.level(ilevel), ...
                    plotdat.level_units);
            end

            myplot(plotdat);

            % create/append to the postscript file
            print(gcf,'-dpsc','-append',psfname{iregion});

            % block to go slow and look at each one ...
            % disp('Pausing, hit any key to continue ...')
            % pause

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

[hprior, prior_legstr] = plot_quantity('prior', plotdat);

nobs_poss       = plotdat.ges_Nposs(:,:,:,plotdat.region);
nobs_used_prior = plotdat.ges_Nused(:,:,:,plotdat.region);
nobs_used_poste = plotdat.anl_Nused(:,:,:,plotdat.region);
anl_Nposs       = sum(plotdat.anl_Nposs(:,:,:,plotdat.region));

if (verbose)
    fprintf('region %d %s level %d nobs_poss %d prior %d poste %d\n', ...
        plotdat.region, plotdat.myvarname, plotdat.mylevel, ...
        sum(nobs_poss), sum(nobs_used_prior), sum(nobs_used_poste))
    fprintf('region %d %s level %d prior_legstr %s\n\n', ...
        plotdat.region, plotdat.myvarname, plotdat.mylevel, prior_legstr)
end

if( isfinite(anl_Nposs) )
    [hposte, poste_legstr] = plot_quantity('posterior', plotdat);
    h = legend([hprior, hposte], prior_legstr, poste_legstr);
else
    h = legend(hprior,prior_legstr);
end

set(h,'Interpreter','none','Box','off','FontSize',figuredata.fontsize)

% Attempt to make plotting robust in the face of 'empty' bins.
% The 't' variable has all the temporal bins specified, 
% so we use that to determine the X axis limits. 

t = plotdat.bincenters;
axlims = [min(t(:)) max(t(:)) plotdat.Yrange];
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

% hokey effort to decide to plot months/days vs. daynum vs.
ttot = plotdat.bincenters(plotdat.Nbins) - plotdat.bincenters(1) + 1;

if ((plotdat.bincenters(1) > 1000) && (ttot > 5))
    datetick('x',6,'keeplimits','keepticks');
    monstr = datestr(plotdat.bincenters(1),21);
    xlabelstring = sprintf('month/day - %s start',monstr);
elseif (plotdat.bincenters(1) > 1000)
    datetick('x',15,'keeplimits')
    monstr = datestr(plotdat.bincenters(1),21);
    xlabelstring = sprintf('%s start',monstr);
else
    xlabelstring = 'days';
end
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

ax2h1 = line(t,nobs_poss,      'Color',figuredata.obs_color,'Parent',ax2);
ax2h2 = line(t,nobs_used_prior,'Color',figuredata.ges_color,'Parent',ax2);
set(ax2h1,'LineStyle','none','Marker','o');
set(ax2h2,'LineStyle','none','Marker',figuredata.ges_marker);

if (sum(nobs_used_poste) > 0)
    ax2h3 = line(t,nobs_used_poste,'Color',figuredata.anl_color,'Parent',ax2);
    set(ax2h3,'LineStyle','none','Marker',figuredata.anl_marker);
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

nevaluated = sum(plotdat.ges_Nqc1(:) + plotdat.ges_Nqc3(:));
if (nevaluated > 0)
    set(get(ax2,'Ylabel'), ...
        'String','# of obs : o=possible, \ast,\diamondsuit=evaluated', ...
        'FontSize', figuredata.fontsize)
else
    set(get(ax2,'Ylabel'), ...
        'String','# of obs : o=possible, \ast,\diamondsuit=assimilated', ...
        'FontSize', figuredata.fontsize)
end


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
if ( ~(isfield(x,'allvarnames') && isfield(x,'allvardims')))
    error('Doh! no ''allvarnames'' and ''allvardims'' components')
end

j = 0;

for i = 1:length(x.allvarnames)
    indx = strfind(x.allvardims{i},'time');
    if (indx > 0)
        j = j + 1;

        basenames{j} = ReturnBase(x.allvarnames{i});
        basedims{ j} = x.allvardims{i};
    end
end

[~,i,~] = unique(basenames);
y     = cell(length(i),1);
ydims = cell(length(i),1);
for k = 1:length(i)
    fprintf('%2d is %s\n',k,basenames{i(k)})
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


function figdata = setfigure()
%%
%  figure out a page layout
%  extra space at the bottom for the date/file annotation
%  extra space at the top because the titles have multiple lines

orientation   = 'landscape';
position      = [0.10 0.15 0.8 0.7];
fontsize      = 16;
linewidth     = 2.5;
obs_color     = [215/255  10/255  83/255]; % obs_red
ges_color     = [  0/255 128/255   0/255]; % prior_green
anl_color     = [  0/255   0/255 255/255]; % poste_blue
rmse_color    = [  0/255   0/255   0/255]; % black
copy_color    = [  0/255 128/255 128/255]; % teal
ges_marker    = '*';
anl_marker    = 'd';
ges_linestyle = '-';
anl_linestyle = '-';

figdata = struct( ...
    'expcolors',  {{'k','r','b','m','g','c','y'}}, ...
    'expsymbols', {{'o','s','d','p','h','s','*'}}, ...
    'prpolines',  {{'-','--'}}, ...
    'position'     , position, ...
    'fontsize'     , fontsize, ...
    'orientation'  , orientation, ...
    'linewidth'    , linewidth, ...
    'obs_color'    , obs_color, ...
    'ges_color'    , ges_color, ...
    'anl_color'    , anl_color, ...
    'rmse_color'   , rmse_color, ...
    'copy_color'   , copy_color, ...
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
    value = ncread(fname,varname);
else
    value = [];
end

%=====================================================================

function [h, legstr] = plot_quantity(phase, plotdat)

global figuredata

switch lower(phase)
    case 'prior'
        data      = plotdat.ges_copy( :,:,:,plotdat.region);
        Nused     = plotdat.ges_Nused(:,:,:,plotdat.region);
        color     = figuredata.ges_color;
        marker    = figuredata.ges_marker;
        linestyle = figuredata.ges_linestyle;
        linewidth = figuredata.linewidth;
        string1   = 'forecast:';
    case 'posterior'
        data      = plotdat.anl_copy( :,:,:,plotdat.region);
        Nused     = plotdat.anl_Nused(:,:,:,plotdat.region);
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
    legstr = sprintf('%s mean = %.5g', string1, data_mean);
else
    legstr = ' ';
end

h = line(plotdat.bincenters,data);
set(h, 'LineStyle',  linestyle, ...
    'LineWidth',  linewidth, ...
    'Color',      color, ...
    'Marker',     marker,    ...
    'MarkerSize', 2*linewidth);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
